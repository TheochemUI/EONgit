#include "SafeNudgedElasticBand.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Log.h"

using namespace helper_functions;

class NEBObjectiveFunction : public ObjectiveFunction
{
    public:

        NEBObjectiveFunction(SafeNudgedElasticBand *nebPassed, Parameters *parametersPassed)
        {
            neb = nebPassed;
            parameters = parametersPassed;
        }

        ~NEBObjectiveFunction(void){};

        VectorXd getGradient(bool fdstep=false)
        {
            VectorXd forceV;
            forceV.resize(3*neb->natoms*neb->nimages);
            if(neb->movedAfterForceCall) neb->updateForces();
            for(size_t counter{0};auto image : neb->imageArray){
                forceV.segment(3*neb->natoms*(counter-1), 3*neb->natoms) =
                    VectorXd::Map(neb->projectedForceArray.at(counter).data(), 3*neb->natoms);
                counter++;
            }
            return -forceV;
        }

        double getEnergy()
        {
            double Energy=0;
            for(long i=1; i<=neb->nimages; i++) {
                Energy += neb->imageArray[i].getPotentialEnergy();
            }
            return Energy;
        }

        void setPositions(VectorXd x)
        {
            neb->movedAfterForceCall = true;
            for(long i=1; i<=neb->nimages; i++) {
                neb->imageArray[i].setPositions(MatrixXd::Map(x.segment(3*neb->natoms*(i-1),3*neb->natoms).data(),neb->natoms,3));
            }
        }

        VectorXd getPositions()
        {
            VectorXd posV;
            posV.resize(3*neb->natoms*neb->nimages);
            for(long i=1; i<=neb->nimages; i++){
                posV.segment(3*neb->natoms*(i-1),3*neb->natoms) = VectorXd::Map(neb->imageArray.at(i).getPositions().data(), 3*neb->natoms);
            }
            return posV;
        }

        int degreesOfFreedom() { return 3*neb->nimages*neb->natoms; }

        bool isConverged() { return getConvergence() < parameters->nebConvergedForce; }

        double getConvergence() { return neb->convergenceForce(); }

        VectorXd difference(VectorXd a, VectorXd b) {
            VectorXd pbcDiff(3*neb->nimages*neb->natoms);
            for (int i=1;i<=neb->nimages;i++) {
                int n = (i-1)*3*neb->natoms;
                int m = 3*neb->natoms;
                pbcDiff.segment(n,m) = neb->imageArray.at(i).pbcV(a.segment(n,m)-b.segment(n,m));
            }
            return pbcDiff;
        }

    private:

        SafeNudgedElasticBand *neb;
        Parameters *parameters;
};


SafeNudgedElasticBand::SafeNudgedElasticBand(Matter initialPassed, Matter finalPassed, Parameters* parametersPassed)
{
    parameters = parametersPassed;
    nimages = parameters->nebImages;
    totImages = nimages + 2;
    natoms = initialPassed.numberOfAtoms();
    log("\nNEB: initialize\n");
    for (size_t idx{0}; idx < totImages; idx++){
        imageArray.emplace_back(initialPassed);
    }
    imageArray.back() = finalPassed;
    tangentArray.resize(totImages);
    projectedForceArray.resize(totImages);
    extremumEnergies.resize(2*(nimages+1));
    extremumPositions.resize(2*(nimages+1));
    extremumCurvatures.resize(2*(nimages+1));
    numExtrema = 0;

    auto posInit = imageArray.front().getPositionsFree();
    auto posFinal = imageArray.back().getPositionsFree();
    const AtomMatrix imageSep = imageArray.front().pbc(posFinal-posInit)/(imageArray.size());
    for (double idx{0}; auto &image: imageArray){
      if (idx == imageArray.size()-1 or (idx == 0)) { // Don't change the final and first image
        ++idx;
        image.setPotential(initialPassed.getPotential());
        continue;
      }
      image.setPositions(posInit + imageSep * idx);
      image.setPotential(initialPassed.getPotential());
      image.getPotentialEnergy();
      // image.getForces();
      // image.useCache = true;
      ++idx;
    }

    movedAfterForceCall = true;
    climbingImage = 0;
    return;
}

SafeNudgedElasticBand::~SafeNudgedElasticBand()
{
    clean();
    return;
}

void SafeNudgedElasticBand::clean(void)
{
}

int SafeNudgedElasticBand::compute(void)
{
    int status = 0;
    long iteration = 0;

    log("Nudged elastic band calculation started.\n");

    updateForces();

    NEBObjectiveFunction objf(this, parameters);

    Optimizer *optimizer = Optimizer::getOptimizer(&objf, parameters);

    const char *forceLabel = parameters->optConvergenceMetricLabel.c_str();
    log("%10s %12s %14s %11s %12s\n", "iteration", "step size", forceLabel, "max image", "max energy");
    log("---------------------------------------------------------------\n");

    char fmt[] = "%10li %12.4e %14.4e %11li %12.4f\n";
    char fmtTiny[] = "%10li %12.4e %14.4e %11li %12.4e\n";

    while (!objf.isConverged())
    {
        if (parameters->writeMovies) {
            bool append = true;
            if (iteration == 0) append = false;
            this->imageArray.at(maxEnergyImage).matter2con("neb_maximage.con", append);
        }
        VectorXd pos = objf.getPositions();
        if(iteration) { // so that we print forces before taking an optimizer step
            if (iteration >= parameters->nebMaxIterations) {
                status = STATUS_BAD_MAX_ITERATIONS;
                break;
            }
            optimizer->step(parameters->optMaxMove);
        }
        iteration++;

        double dE = this->imageArray[maxEnergyImage].getPotentialEnergy() -
                    imageArray.front().getPotentialEnergy();
        std::cout<<dE<<" is the energy diff\n";
        double stepSize = helper_functions::maxAtomMotionV(imageArray.front().pbcV(objf.getPositions()-pos));
        if (dE > 0.01) {
            log(fmt, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
        }else{
            log(fmtTiny, iteration, stepSize, convergenceForce(), maxEnergyImage, dE);
        }
    }

    if(objf.isConverged()) {
        status = STATUS_GOOD;
        log("NEB converged\n");
    }

    printImageData();
    findExtrema();

    delete optimizer;
    return status;
}


// generate the force value that is compared to the convergence criterion
double SafeNudgedElasticBand::convergenceForce(void)
{
    if(movedAfterForceCall) updateForces();
    double fmax = 0;

    for(long i=1; i<=nimages; i++) {

            if( parameters->nebClimbingImageConvergedOnly == true &&
                parameters->nebClimbingImageMethod &&
                climbingImage!=0 ) {
                i = climbingImage;
            }
            if (parameters->optConvergenceMetric == "norm") {
                fmax = max(fmax, projectedForceArray[i].norm());
            } else if (parameters->optConvergenceMetric == "max_atom") {
                for (int j=0;j<imageArray.front().numberOfAtoms();j++) {
                    if (imageArray.front().getFixed(j)) continue;
                    fmax = max(fmax, projectedForceArray[i].row(j).norm());
                }
            } else if (parameters->optConvergenceMetric == "max_component") {
                fmax = max(fmax, projectedForceArray[i].maxCoeff());
            } else {
                log("[Nudged Elastic Band] unknown opt_convergence_metric: %s\n", parameters->optConvergenceMetric.c_str());
                exit(1);
            }
            if( parameters->nebClimbingImageConvergedOnly == true &&
                parameters->nebClimbingImageMethod &&
                climbingImage!=0 ) {
                break;
            }
    }
    return fmax;
}

// Update the forces, do the projections, and add spring forces
void SafeNudgedElasticBand::updateForces(void)
{
    // variables for tangent
    double maxDiffEnergy, minDiffEnergy;
    double energyDiffPrev, energyDiffNext;
    double energy, energyPrev, energyNext;
    //bool higherEnergyPrev, higherEnergyNext;

    // variables for climbing image
    double maxEnergy;

    // variables for force projections
    AtomMatrix force(natoms,3), forcePerp(natoms,3), forcePar(natoms,3);
    AtomMatrix forceSpringPar(natoms,3), forceSpring(natoms,3), forceSpringPerp(natoms,3);
    AtomMatrix forceDNEB(natoms, 3);
    AtomMatrix pos(natoms,3), posNext(natoms,3), posPrev(natoms,3);
    double distNext, distPrev;

    // update the forces on the images and find the highest energy image
    maxEnergy = this->imageArray.front().getPotentialEnergy();
    std::cout<<"max energy is "<<maxEnergy<<"\n";
    maxEnergyImage = 0;
    for (size_t idx{1}; auto &image: imageArray){
      if (idx == imageArray.size()-1) {
        ++idx;
        continue;
      }
      image.getForces();
      double maybeMaxE = image.getPotentialEnergy();
        if(maybeMaxE > maxEnergy) {
            maxEnergy = maybeMaxE;
            std::cout<<"max energy is now "<<maxEnergy<<"\n";
            maxEnergyImage = idx;
        }
      ++idx;
    }

    for (size_t idx{1}; auto &image: imageArray){
      if (idx == imageArray.size()-1) {
        ++idx;
        continue;
      }
        // set local variables
        force = image.getForces();
        pos = image.getPositions();
        posPrev = imageArray.at(idx-1).getPositions();
        posNext = imageArray.at(idx+1).getPositions();
        std::cout<<"Within loop "<<idx<<" has "<<imageArray.at(idx).getPotentialEnergy()<<" energy\n";
        energy = imageArray[idx].getPotentialEnergy();
        energyPrev = imageArray.at(idx-1).getPotentialEnergy();
        energyNext = imageArray.at(idx+1).getPotentialEnergy();

        // determine the tangent
        if(parameters->nebOldTangent) {
            // old tangent
            tangentArray.at(idx) = image.pbc(posNext - posPrev);
        }else{
          // new improved tangent
	  //higherEnergyPrev = energyPrev > energyNext;
	  //higherEnergyNext = energyNext > energyPrev;

	  if(energyNext > energy && energy > energyPrev) {
	    tangentArray.at(idx) = image.pbc(posNext - pos);
	  }else if(energy > energyNext && energyPrev > energy){
	    tangentArray.at(idx) = image.pbc(pos - posPrev);
	  }else{
	    // we are at an extremum
	    energyDiffPrev = energyPrev - energy;
	    energyDiffNext = energyNext - energy;

	    // calculate the energy difference to neighboring images
	    minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
	    maxDiffEnergy = max(abs(energyDiffPrev), abs(energyDiffNext));

	    // use these energy differences to weight the tangent
	    if(energyDiffPrev > energyDiffNext) {
	      tangentArray.at(idx) = image.pbc(posNext - pos) * minDiffEnergy;
	      tangentArray.at(idx) += image.pbc(pos - posPrev) * maxDiffEnergy;
	    }else{
	      tangentArray.at(idx) = image.pbc(posNext - pos) * maxDiffEnergy;
	      tangentArray.at(idx) += image.pbc(pos - posPrev) * minDiffEnergy;
	    }
	  }
	}
        tangentArray.at(idx).normalize();

        // project the forces and add springs
        force = image.getForces();

        // calculate the force perpendicular to the tangent
        forcePerp = force - (force.array() * (tangentArray.at(idx)).array()).sum() * tangentArray.at(idx);
        forceSpring = parameters->nebSpring * image.pbc((posNext - pos) - (pos - posPrev));

        // calculate the spring force
        distPrev = image.pbc(posPrev - pos).squaredNorm();
        distNext = image.pbc(posNext - pos).squaredNorm();
        forceSpringPar = parameters->nebSpring * (distNext-distPrev) * tangentArray.at(idx);

        if (parameters->nebDoublyNudged) {
            forceSpringPerp = forceSpring - (forceSpring.array() * (tangentArray.at(idx)).array()).sum() * tangentArray.at(idx);
            forceDNEB = forceSpringPerp - (forceSpringPerp.array() * forcePerp.normalized().array()).sum() * forcePerp.normalized();
            if (parameters->nebDoublyNudgedSwitching) {
                double switching;
                switching = 2.0/M_PI * atan(pow(forcePerp.norm(),2) / pow(forceSpringPerp.norm(),2));
                forceDNEB *= switching;
            }
        }else{
            forceDNEB.setZero();
        }


        if(parameters->nebClimbingImageMethod && idx==maxEnergyImage)
        {
            // we are at the climbing image
            climbingImage = maxEnergyImage;
            projectedForceArray.at(idx) = force - (2.0 * (force.array() * (tangentArray.at(idx)).array()).sum() * tangentArray.at(idx)) + forceDNEB;
        }
        else  // all non-climbing images
        {
            // sum the spring and potential forces for the neb force
            if (parameters->nebElasticBand) {
                projectedForceArray.at(idx) = forceSpring + force;
            }else{
                projectedForceArray.at(idx) = forceSpringPar + forcePerp + forceDNEB;
            }
            //projectedForceArray.at(idx) = forceSpring + forcePerp;


            //if (parameters->nebFullSpring) {


            movedAfterForceCall = false;  // so that we don't repeat a force call
        }

        //zero net translational force
        if (imageArray[idx].numberOfFreeAtoms() == imageArray[idx].numberOfAtoms()) {
            for (int j=0;j<=2;j++) {
                double translationMag = projectedForceArray.at(idx).col(j).sum();
                int natoms = projectedForceArray.at(idx).col(j).size();
                projectedForceArray.at(idx).col(j).array() -= translationMag/((double)natoms);
            }
        }
    }

    return;
}

// Print NEB image data
void SafeNudgedElasticBand::printImageData(bool writeToFile)
{
    double dist, distTotal=0;
    AtomMatrix tangentStart = imageArray[0].pbc(imageArray[1].getPositions() - imageArray[0].getPositions());
    AtomMatrix tangentEnd = imageArray[nimages].pbc(imageArray[nimages+1].getPositions() - imageArray[nimages].getPositions());
    AtomMatrix tang;

    log("Image data (as in neb.dat)\n");

    FILE *fh=NULL;
    if (writeToFile) {
        fh = fopen("neb.dat", "w");
    }

    for(long i=0; i<=nimages+1; i++)
    {
        if(i==0){ tang = tangentStart; }
        else if (i==nimages+1) { tang = tangentEnd; }
        else { tang = tangentArray[i]; }
        if(i>0) {
            dist = imageArray[i].distanceTo(imageArray[i-1]);
            distTotal += dist;
        }
        if (fh == NULL) {
            log("%3li %12.6f %12.6f %12.6f\n",i,distTotal,
                imageArray[i].getPotentialEnergy()-imageArray[0].getPotentialEnergy(), (imageArray[i].getForces().array()*tang.array()).sum());
        }else{
            fprintf(fh, "%3li %12.6f %12.6f %12.6f\n",i,distTotal,
                imageArray[i].getPotentialEnergy()-imageArray[0].getPotentialEnergy(), (imageArray[i].getForces().array()*tang.array()).sum());
        }
    }
    if (writeToFile) {
        fclose(fh);
    }
}

// Estimate the barrier using a cubic spline
void SafeNudgedElasticBand::findExtrema(void)
{
    // calculate the cubic parameters for each interval (a,b,c,d)

    AtomMatrix tangentEndpoint;
    double a[nimages+1], b[nimages+1], c[nimages+1], d[nimages+1];
    double F1, F2, U1, U2, dist;

    for(long i=0; i<=nimages; i++)
    {
        dist = imageArray[i].distanceTo(imageArray[i+1]);
        if(i==0) {
            tangentEndpoint = imageArray[i].pbc(imageArray[1].getPositions() - imageArray[0].getPositions());
            tangentEndpoint.normalize();
            F1 = (imageArray[i].getForces().array()*tangentEndpoint.array()).sum()*dist;
        } else {
            F1 = (imageArray[i].getForces().array()*(tangentArray[i]).array()).sum()*dist;
        }
        if(i==nimages) {
            tangentEndpoint =  imageArray[i+1].pbc(imageArray[nimages+1].getPositions() - imageArray[nimages].getPositions());
            tangentEndpoint.normalize();
            F2 = (imageArray[i+1].getForces().array()*tangentEndpoint.array()).sum()*dist;
        } else {
            F2 = (imageArray[i+1].getForces().array()*(tangentArray[i+1]).array()).sum()*dist;
        }
        U1 = imageArray[i].getPotentialEnergy();
        U2 = imageArray[i+1].getPotentialEnergy();
        a[i] = U1;
        b[i] = -F1;
        c[i] = 3.*(U2 - U1) + 2.*F1 + F2;
        d[i] = -2.*(U2 - U1) - (F1 + F2);
    }

    // finding extrema along the MEP

//    long numExtrema = 0;
//    double extremaEnergy[2*(images+1)]; // the maximum number of extrema
//    double extremaPosition[2*(images+1)];
    double discriminant, f;

    for(long i=0; i<=nimages; i++)
    {
        discriminant = pow(c[i],2) - 3.*b[i]*d[i];
        if(discriminant >= 0) {
            f = -1;

            // quadratic case
            if( (d[i] == 0) && (c[i] != 0) ) {
                f = ( -b[i]/(2.*c[i]) );
            }
            // cubic case 1
            else if( d[i] != 0 ) {
                f = -(c[i] + sqrt(discriminant))/(3.*d[i]);
            }
            if( (f >= 0) && (f <= 1) ) {
                extremumPositions[numExtrema] = i + f;
                extremumEnergies[numExtrema] = d[i]*pow(f,3) + c[i]*pow(f,2) + b[i]*f + a[i];
                extremumCurvatures[numExtrema] = 6.0*d[i]*f + 2*c[i];
                numExtrema ++;
            }
            // cubic case 2
            if( d[i] != 0 ) {
                f = ( -(c[i] - sqrt(discriminant))/(3.*d[i]) );
            }
            if( (f >= 0) && (f <= 1) ) {
                extremumPositions[numExtrema] = i + f;
                extremumEnergies[numExtrema] = d[i]*pow(f,3) + c[i]*pow(f,2) + b[i]*f + a[i];
                extremumCurvatures[numExtrema] = 6*d[i]*f + 2*c[i];
                numExtrema ++;
            }
        }
    }

    log("\nFound %li extrema\n",numExtrema);
    log("Energy reference: %f\n",imageArray.front().getPotentialEnergy());
    for(long i=0; i<numExtrema; i++) {
        log("extrema #%li at image position %f with energy %f and curvature %f\n",i+1,extremumPositions[i],extremumEnergies[i]-imageArray.front().getPotentialEnergy(), extremumCurvatures[i]);
    }
}
