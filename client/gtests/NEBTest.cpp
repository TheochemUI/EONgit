/*
 * Test.cpp
 *
 *  Created on: 23 Feb 2022
 *      Author: Rohit Goswami
 *     Company: University of Iceland
 */

#include "NEBTest.h"
#include "../potentials/Morse/Morse.cpp"
#include "Eigen/src/Core/Matrix.h"
#include "NudgedElasticBand.h"

#include <algorithm>

using namespace std::placeholders;

namespace tests {

NEBTest::NEBTest()
    : params{new Parameters},
      m1{new Matter(params)},
      m2{new Matter(params)},
      threshold{1e-3} {}

NEBTest::~NEBTest() {
  delete params;
  delete m1;
}

void NEBTest::SetUp() {
  std::string confile("pos.con");
  m1->con2matter(confile);
  m2->con2matter(confile);
  m2->setPositions(m1->getPositions() * 0.7);
  params->optMaxIterations = 100000;
}

void NEBTest::TearDown() {}

TEST_F(NEBTest, TestCreation) {
  NudgedElasticBand *neb = new NudgedElasticBand(m1, m2, params);
  NudgedElasticBand::NEBStatus status = neb->compute();
  ASSERT_EQ(status, NudgedElasticBand::NEBStatus::STATUS_GOOD);
  ASSERT_EQ(neb->atoms, m1->numberOfAtoms());
  ASSERT_EQ(neb->images, 5);
  ASSERT_EQ(neb->climbingImage, 0);
  ASSERT_EQ(neb->numExtrema, 4);
  ASSERT_EQ(neb->movedAfterForceCall, false);
  ASSERT_NEAR(*neb->extremumEnergy, -38.540506834483423, threshold);
  ASSERT_NEAR(*neb->extremumPosition, 2.7250638348703475, threshold);
  ASSERT_NEAR(*neb->extremumCurvature, -11.971787908472798, threshold);
  ASSERT_EQ(neb->maxEnergyImage, 6);
}

TEST_F(NEBTest, TestObjective) {
  auto vecEq =
      std::bind(helper_functions::eigenEquality<VectorType>, _1, _2, 1e-4);
  VectorType getPositions_{
      {9.95383, 11.7605, 12.2282, 9.3464,  11.462,  10.0695, 12.9505, 12.3277,
       11.0109, 11.9457, 12.7707, 10.2873, 13.1479, 11.0061, 11.1346, 12.3646,
       10.1611, 10.5346, 13.4712, 12.9691, 11.6451, 11.6785, 13.7695, 10.3685,
       11.7904, 12.4083, 9.32115, 12.3438, 9.16952, 10.8255, 13.8191, 10.6776,
       11.8606, 11.9621, 10.3164, 9.59155, 10.4243, 11.5083, 10.9473, 9.42994,
       11.1415, 11.5846, 8.85448, 10.8587, 9.5395,  12.2689, 11.6789, 10.4314,
       11.317,  12.0986, 9.74591, 12.4559, 10.4268, 10.5486, 11.7139, 9.62629,
       9.98012, 12.7621, 12.2865, 11.0322, 11.0638, 13.0448, 9.82279, 11.1699,
       11.7553, 8.83056, 11.6941, 8.68691, 10.2558, 13.0917, 10.1156, 11.2364,
       11.3326, 9.77343, 9.08673, 9.87566, 10.9026, 10.3711, 8.90606, 10.5225,
       10.941,  8.36257, 10.2555, 9.00953, 11.5873, 11.0301, 9.85189, 10.6882,
       11.4264, 9.20447, 11.7639, 9.84752, 9.96258, 11.0631, 9.09149, 9.42567,
       12.0531, 11.6039, 10.4193, 10.4492, 12.3201, 9.27708, 10.5493, 11.1022,
       8.33998, 11.0444, 8.2043,  9.68601, 12.3644, 9.5536,  10.6121, 10.703,
       9.23046, 8.58192, 9.32701, 10.2969, 9.79497, 8.38217, 9.90357, 10.2974,
       7.87065, 9.65221, 8.47956, 10.9057, 10.3812, 9.27237, 10.0595, 10.7543,
       8.66303, 11.0719, 9.26826, 9.37654, 10.4123, 8.5567,  8.87122, 11.3441,
       10.9214, 9.80643, 9.83451, 11.5954, 8.73137, 9.92878, 10.4491, 7.84939,
       10.3948, 7.7217,  9.11625, 11.6371, 8.99162, 9.9879,  10.0734, 8.68749,
       8.0771,  8.77836, 9.69117, 9.2188,  7.85828, 9.2846,  9.65383, 7.37874,
       9.04894, 7.94958, 10.2241, 9.73241, 8.69285, 9.4308,  10.0821, 8.12159,
       10.3799, 8.68899, 8.79051, 9.76155, 8.02191, 8.31677, 10.6351, 10.2388,
       9.19352, 9.21985, 10.8707, 8.18566, 9.30823, 9.79606, 7.3588,  9.7451,
       7.23909, 8.54648, 10.9098, 8.42965, 9.36366, 9.4438,  8.14452, 7.57228,
       8.22971, 9.08547, 8.64262}};
  VectorType getGradient_{
      {-0.620969, 0.324815,   1.72152,   -1.43938,  -0.0774983, -1.16484,
       6.6729,    7.67128,    7.38725,   -6.05536,  8.41656,    -9.53704,
       8.4028,    -3.75277,   8.44217,   -6.06126,  -11.6783,   -8.01197,
       -4.91462,  -6.07976,   -6.04805,  2.73872,   -10.2805,   -0.682081,
       1.42192,   3.21242,    8.495,     0.262068,  11.1458,    -3.13862,
       -6.27767,  3.0795,     -6.84755,  4.24387,   -1.64951,   9.76833,
       1.62699,   -0.33197,   -0.384115, -0.818931, 0.409246,   2.19274,
       -1.86926,  -0.102853,  -1.48183,  17.1804,   20.647,     20.0654,
       -14.8618,  22.7995,    -26.4949,  21.8165,   -10.4645,   22.8444,
       -15.1498,  -30.7575,   -21.5905,  -15.0624,  -18.6022,   -18.4567,
       8.30498,   -31.0895,   -2.32524,  4.36909,   10.0838,    26.6817,
       0.739972,  33.0106,    -9.517,    -19.2817,  9.44784,    -20.9579,
       12.7003,   -4.95337,   29.48,     1.93258,   -0.427957,  -0.440113,
       -1.02858,  0.467385,   2.56483,   -2.2771,   -0.130767,  -1.7286,
       42.0194,   51.6176,    50.3819,   -35.4384,  57.1738,    -67.2247,
       53.5858,   -26.5962,   57.2576,   -36.4835,  -76.0878,   -54.008,
       -39.8488,  -49.1812,   -48.7453,  21.8785,   -81.8109,   -6.39807,
       11.5705,   26.9277,    71.2668,   1.89131,   86.1262,    -25.0707,
       -51.0551,  25.0053,    -55.4161,  33.2784,   -12.9975,   77.5111,
       1.9076,    -0.513508,  -0.390689, -1.12592,  0.383669,   2.28825,
       -2.29861,  -0.147486,  -1.52838,  102.083,   126.979,    124.241,
       -84.8123,  140.889,    -166.734,  130.504,   -66.0263,   141.055,
       -87.8316,  -186.065,   -132.903,  -101.051,  -124.674,   -123.503,
       55.3594,   -206.889,   -16.5456,  29.3565,   68.6107,    181.606,
       4.70967,   216.833,    -63.4359,  -129.526,  63.4233,    -140.488,
       83.9704,   -32.8202,   195.931,   0.662049,  -0.496228,  0.0171093,
       -0.681294, -0.203186,  -0.382335, -0.730815, -0.104022,  0.313015,
       252.606,   316.629,    310.256,   -207.902,  351.679,    -417.821,
       323.424,   -165.548,   352.035,   -216.109,  -462.279,   -331.467,
       -256.552,  -316.466,   -313.394,  140.371,   -524.418,   -42.4782,
       74.554,    174.67,     462.367,   11.8302,   548.196,    -160.849,
       -328.93,   161.041,    -356.621,  212.573,   -83.1205,   496.518,
       -4.45161,  -0.0761206, 1.52368}};
  VectorType diffDat_{
      {-10.5748, -11.4357, -10.5067, -10.7858, -11.5395, -11.2343, -6.2776,
       -4.65642, -3.62365, 6.99894,  -4.35414, 5.17566,  -4.7451,  10.2411,
       -2.69243, 6.57414,  3.1606,   6.45343,  6.61418,  5.95114,  7.30685,
       -8.93978, 0.95,     -11.0506, -10.3685, -9.19588, -0.82615, -12.0817,
       1.97628,  11.0359,  4.90323,  -7.5981,  6.29185,  -7.71823, -11.9659,
       0.17678,  -8.79731, -11.8403, -11.3314, -10.2489, -10.7323, -9.39186,
       -10.7237, -10.9616, -11.0213, 4.9115,   8.9681,   9.634,    -1.1788,
       10.7009,  -11.2408, 9.3606,   4.1087,   12.2958,  -1.8637,  9.61621,
       -6.57062, -2.8245,  -5.8887,  -4.4889,  -2.75882, 5.8657,   -12.148,
       -6.80081, -1.6715,  -7.14886, -10.9541, -0.67631, 5.2272,   -7.3734,
       -0.66776, -7.1943,  1.3677,   10.2732,  -4.60673, -7.94308, -11.3306,
       -10.8112, -9.93464, -10.0551, -8.37617, -10.6397, -10.3863, -10.7381,
       5.4321,   -9.4125,  -9.46999, 3.8734,   -4.2526,  -1.42917, -8.1781,
       -11.4437, -2.70498, 2.4534,   -10.1793, 11.5663,  -1.9019,  -10.7851,
       -9.1646,  11.4293,  5.869,    9.32485,  1.0212,   -9.1745,  -12.0732,
       -9.15309, 2.9219,   -9.75671, 11.5805,  -9.5483,  8.9718,   -2.4246,
       2.77204,  -6.07082, -7.41941, -10.8104, -10.1857, -9.50809, -9.5199,
       -8.00915, -10.1693, -9.7997,  -10.0079, -8.8227,  -8.4022,  -10.0314,
       5.1282,   5.1347,   -0.39703, -5.5679,  -0.29456, 6.67846,  1.7561,
       5.3783,   8.22578,  -12.3951, -10.5954, -8.30943, -4.47511, 6.5156,
       -0.27697, -5.57228, 8.1616,   -1.24339, -5.68513, 9.1113,   2.44785,
       8.8369,   4.43168,  -0.4759,  -1.103,   8.49231,  -12.1461, -8.11631,
       -10.1874, -9.20169, -8.53957, -9.48779, -10.0362, -8.10956, -9.15296,
       -7.63657, -7.6181,  6.89659,  1.56315,  7.6672,   -8.4031,  -0.94259,
       -11.9559, 0.76301,  -6.75551, -0.87055, 4.69909,  10.2162,  7.8129,
       -1.7048,  2.41248,  6.15115,  -10.2887, -0.66386, -9.75423, -10.1261,
       5.0082,   2.0851,   -9.04309, 5.60452,  10.1602,  2.61135,  9.01534,
       3.1292,   8.73498,  -11.0543, 12.3187,  -9.16159, -7.11894}};
  VectorType setPositionDat_{VectorType::Zero(getPositions_.size())};
  NudgedElasticBand *neb = new NudgedElasticBand(m1, m2, params);
  NEBObjectiveFunction nebo{neb, params};
  ASSERT_NEAR(nebo.getEnergy(), 231.60596556863382, threshold);
  ASSERT_EQ(nebo.isConverged(), false);
  ASSERT_EQ(nebo.degreesOfFreedom(), 195);
  ASSERT_NEAR(nebo.getConvergence(), 1710.2806025488981, threshold);
  ASSERT_PRED2(vecEq, getPositions_, nebo.getPositions());
  ASSERT_PRED2(vecEq, getGradient_, nebo.getGradient());
  VectorType nebo_diff = nebo.difference(getGradient_, getPositions_);
  ASSERT_PRED2(vecEq, diffDat_, nebo_diff);
  // Test setting
  nebo.setPositions(VectorType::Ones(getPositions_.size()).array() - 1);
  ASSERT_PRED2(vecEq, setPositionDat_, nebo.getPositions());
}

} /* namespace tests */
