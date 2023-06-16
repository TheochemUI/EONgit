#include "CommandLine.h"
#include "version.h"

void singlePoint(std::unique_ptr<Matter> matter) {
  fmt::printf("Energy:         %.10f\n", matter->getPotentialEnergy());
  fmt::print("(free) Forces:         \n{}\n", fmt::streamed(matter->getForcesFree()));
  fmt::printf("Max atom force: %.10g\n", matter->maxForce());
}

void minimize(std::unique_ptr<Matter> matter, string confileout) {
  matter->relax(false, false);
  if (confileout.length() > 0) {
    fmt::printf("saving relaxed structure to %s\n", confileout.c_str());
  } else {
    fmt::printf("no output file specified, not saving\n");
  }
  matter->matter2con(confileout);
}

void usage(void) {
  fmt::fprintf(stderr, "Usage: eonclient [options] inputConfile [outputConfile]\n");
  char fmtStr[] = "  -%-2s %s\n";

  fmt::fprintf(stderr, "Job Type:\n");
  fmt::fprintf(stderr, fmtStr, "v", "print version information");
  fmt::fprintf(stderr, fmtStr, "m",
          "Minimization of inputConfile saves to outputConfile");
  fmt::fprintf(stderr, fmtStr, "s", "Single point energy of inputConfile");
  fmt::fprintf(stderr, fmtStr, "c",
          "Compare structures of inputConfile to outputConfile");
  fmt::fprintf(stderr, fmtStr, "o", "Optimization method [default: qm]");
  fmt::fprintf(stderr, fmtStr, "f", "Convergence force [default: 0.001]");
  fmt::fprintf(stderr, fmtStr, "t", "Distance tolerance [default: 0.1]");

  fmt::fprintf(stderr, "Required Options:\n");
  fmt::fprintf(stderr, fmtStr, "p", "The potential (e.g. qsc, lj, eam_al)");
}

void commandLine(int argc, char **argv) {
// no getopt on windows
#ifndef WIN32
  int c;
  bool sflag = false, mflag = false, pflag = false, cflag = false;
  double optConvergedForce = 0.001;

  string potential;
  string confile;
  string optimizer("cg");

  auto params = std::make_shared<Parameters>();

  while ((c = getopt(argc, argv, "chsmp:f:o:t:v")) != -1) {
    switch (c) {
    case 'c':
      cflag = true;
      break;
    case 's':
      sflag = true;
      break;
    case 'm':
      mflag = true;
      break;
    case 'p':
      pflag = true;
      potential = optarg;
      break;
    case 't':
      params->distanceDifference = atof(optarg);
      break;
    case 'o':
      optimizer = optarg;
      break;
    case 'f':
      cout << optarg << endl;
      optConvergedForce = atof(optarg);
      break;
    case 'h':
      usage();
      exit(0);
    case 'v':
      fmt::printf("eonclient version r%s\n", VERSION);
      fmt::printf("          compiled %s\n", BUILD_DATE);
      exit(0);
    case '?':
      if (optopt == 'p')
        fmt::fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else
        fmt::fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      usage();
      exit(2);
    }
  }

  if (sflag && mflag) {
    fmt::fprintf(stderr, "Cannot specify both minimization and single point\n");
    exit(2);
  }

  if (!pflag && (sflag || mflag)) {
    fmt::fprintf(stderr, "Must specify a potential\n");
    exit(2);
  } else if (!cflag) {
    for (string::size_type i = 0; i < potential.length(); ++i) {
      potential[i] = tolower(potential[i]);
    }
  }

  int extraArgs = argc - optind;

  if (extraArgs < 1) {
    fmt::fprintf(stderr, "Only one non-option argument is allowed: the con file\n");
    exit(2);
  } else {
    confile = argv[optind];
  }

  if (!cflag) {
    params->potential = helper_functions::getPotentialType(potential);
  }
  params->optMethod = optimizer;
  params->optConvergedForce = optConvergedForce;

  // Logger::getInstance().logInit("log.txt");

  auto pot = helper_functions::makePotential(params);
  auto matter = std::make_unique<Matter>(pot, params);
  auto matter2 = std::make_unique<Matter>(pot, params);
  matter->con2matter(confile);

  string confileout;
  if (extraArgs == 2) {
    confileout = argv[optind + 1];
    if (cflag)
      matter2->con2matter(confileout);
  }

  if (sflag) {
    singlePoint(std::move(matter));
  } else if (mflag) {
    minimize(std::move(matter), confileout);
  } else if (cflag) {
    params->checkRotation = true;
    if (matter->compare(*matter2, true)) {
      fmt::printf("structures match\n");
    } else {
      fmt::printf("structures do not match\n");
    }
  }

#endif
}
