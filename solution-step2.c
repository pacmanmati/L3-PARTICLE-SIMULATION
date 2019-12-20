// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>

double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;
int NumberOfBodies = 0;
double** x; // Pointer to pointers. Each pointer in turn points to three coordinates, i.e. each pointer represents one molecule/particle/body.
double** v; //  Equivalent to x storing the velocities.
double*  mass; // One mass entry per molecule/particle.
double   timeStepSize = 0.0; // Global time step size used.
double   maxV; // Maximum velocity of all particles.
double   minDx; // Minimum distance between two elements.

/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}

std::ofstream videoFile;
/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}

/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }
  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}

/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  double* distances = new double[(NumberOfBodies * (NumberOfBodies - 1)) / 2];

  double* force0 = new double[NumberOfBodies]; // force along x direction
  double* force1 = new double[NumberOfBodies]; // force along y direction
  double* force2 = new double[NumberOfBodies]; // force along z direction

  double diameter = 0.2;

  // check if a collision occured
  for (int i = 0; i < NumberOfBodies;  i++) {
    for (int j = 0; j < NumberOfBodies; j++) {
      if (i == j) continue; // skip same comparison
      /* std::cout << "distance on x: " << std::abs(x[i][0] - x[j][0]) << "diameter: " << diameter << std::endl; */
      // check if bodies will collide
      if (std::abs(x[i][0] - x[j][0]) < diameter &&
	  std::abs(x[i][1] - x[j][1]) < diameter &&
	  std::abs(x[i][2] - x[j][2]) < diameter) {
	/* std::cout << "ja ja collision" << std::endl; */
	int first, last;
	if (i < j) {
	  first = i;
	  last = j;
	} else {
	  first = j;
	  last = i;
	}
	// update velocity and position of the first index particle
	v[first][0] = v[i][0]*mass[i]/(mass[i]+mass[j]) +
	  v[j][0]*mass[j]/(mass[i]+mass[j]);
	v[first][1] = v[i][1]*mass[i]/(mass[i]+mass[j]) +
	  v[j][1]*mass[j]/(mass[i]+mass[j]);
	v[first][2] = v[i][2]*mass[i]/(mass[i]+mass[j]) +
	  v[j][2]*mass[j]/(mass[i]+mass[j]);
	
	x[first][0]= x[i][0]*mass[i]/(mass[i]+mass[j]) +
	  x[j][0]*mass[j]/(mass[i]+mass[j]);
	x[first][1]= x[i][1]*mass[i]/(mass[i]+mass[j]) +
	  x[j][1]*mass[j]/(mass[i]+mass[j]);
	x[first][2] = x[i][2]*mass[i]/(mass[i]+mass[j]) +
	  x[j][2]*mass[j]/(mass[i]+mass[j]);
	
	mass[first] += mass[last];
	// shift values s.t. positions, masses and velocities are maintained, allowing for safe removal of the redundant particle
	for (int k = last+1; k < NumberOfBodies; k++) {
	  v[k-1] = v[k];
	  x[k-1][0] = x[k][0];
	  x[k-1][1] = x[k][1];
	  x[k-1][2] = x[k][2];
	  mass[k-1] = mass[k];
	}
	NumberOfBodies--;
      }
    }
  }

  // if we're on the last particle
  if (NumberOfBodies == 1) {
    std::cout << "last body @ (X, Y, Z) : (" << x[0][0] << ", " << x[0][1] << ", "<< x[0][2] << ")" << std::endl;
    std::exit(0);
  }

  // zero all the forces every iteration
  for (int i = 0; i < NumberOfBodies; i++) {
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;
  }

  int distStore = 0;
  int distFetch = 0;

  for (int j = 0; j < NumberOfBodies; j++) { // work out the jth particles forces
    for (int i = 0; i < NumberOfBodies; i++) { // iterate over (i-1) other particles
      if (i == j) continue; // skip calculating particle's force on itself
      double distance;
      // if i < j then we have already cached the result
      if (i < j) {
	distance = distances[distFetch];
	distFetch++;
      } else {
	// uncached: find distance between the jth and ith particle
	distance = sqrt((x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
				     (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
				     (x[j][2]-x[i][2]) * (x[j][2]-x[i][2]));
	distances[distStore] = distance;
	distStore++;
      }
      // x,y,z forces acting on particle j
      force0[j] += (x[i][0]-x[j][0]) * mass[i]*mass[j] / distance / distance / distance;
      force1[j] += (x[i][1]-x[j][1]) * mass[i]*mass[j] / distance / distance / distance;
      force2[j] += (x[i][2]-x[j][2]) * mass[i]*mass[j] / distance / distance / distance;
      // save minDx
      minDx = std::min(minDx, distance);
    }
  }

  // update position of every particle using velocity
  for (int i = 0; i < NumberOfBodies; i++) {
    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];

    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

    maxV = std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  }

  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      /* std::cout << "plot next snapshot" */
      /* 		    << ",\t time step=" << timeStepCounter */
      /* 		    << ",\t t="         << t */
      /* 				<< ",\t dt="        << timeStepSize */
      /* 				<< ",\t v_max="     << maxV */
      /* 				<< ",\t dx_min="    << minDx */
      /* 				<< std::endl; */

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
