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
  /* double **forces = new double*[NumberOfBodies]; // cleanup forces into one pointer for nicer code */
  /* double *fpool = new double[3*NumberOfBodies](); */
  /* for (int i = 0; i < NumberOfBodies; i++, fpool+=3) forces[i] = fpool; */
  /* auto forces = new double[3][n](); */
  // no idea how to make the single forces pointer work with parallelisation for this step so it's a bit different from step 1
  double *force0 = new double[NumberOfBodies]();
  double *force1 = new double[NumberOfBodies]();
  double *force2 = new double[NumberOfBodies]();
  // appending empty () to the new operator initialises with 0
  /* for (int i = 0; i < NumberOfBodies; i++) forces[i] = new double[3](); */
  // convention: always arr[i][dim]

#pragma omp parallel for reduction(+:force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies]) reduction(min:minDx) schedule(dynamic, 1)
  for (int j = 0; j < NumberOfBodies; j++) { // work out the jth particles forces
/* #pragma omp parallel for reduction(+:forces[j:NumberOfBodies][:3]) reduction(min:minDx) schedule(dynamic, 1) */
    for (int i = j+1; i < NumberOfBodies; i++) { // iterate over (i-1) other particles
      double distanceSquared = 0;
      // find distance between the jth and ith particle
      distanceSquared += (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) + (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) + (x[j][2]-x[i][2]) * (x[j][2]-x[i][2]);
      /* for (int dim = 0; dim < 3; dim++) distanceSquared += (x[j][dim]-x[i][dim]) * (x[j][dim]-x[i][dim]); */
      double distance = std::sqrt(distanceSquared);
      force0[i] += (x[j][0]-x[i][0]) * mass[i]*mass[j] / distance / distance / distance;
      force1[i] += (x[j][1]-x[i][1]) * mass[i]*mass[j] / distance / distance / distance;
      force2[i] += (x[j][2]-x[i][2]) * mass[i]*mass[j] / distance / distance / distance;
      force0[j] -= (x[j][0]-x[i][0]) * mass[i]*mass[j] / distance / distance / distance;
      force1[j] -= (x[j][1]-x[i][1]) * mass[i]*mass[j] / distance / distance / distance;
      force2[j] -= (x[j][2]-x[i][2]) * mass[i]*mass[j] / distance / distance / distance;
      // save minDx
      minDx = std::min(minDx, distance);
    }
  }
  // once we have all of our forces:
  // update pos and vel of each particle
#pragma omp parallel for reduction(max:maxV) schedule(dynamic, 1)
  for (int i = 0; i < NumberOfBodies; i++) {
    x[i][0] += timeStepSize *v[i][0];
    x[i][1] += timeStepSize *v[i][1];
    x[i][2] += timeStepSize *v[i][2];

    v[i][0] += timeStepSize * force0[i] / mass[i];
    v[i][1] += timeStepSize * force1[i] / mass[i];
    v[i][2] += timeStepSize * force2[i] / mass[i];
    /* for (int dim = 0; dim < 3; dim++) { */
    /*   x[i][dim] += timeStepSize * v[i][dim]; */
    /*   v[i][dim] += timeStepSize * forces[i][dim] / mass[i]; */
    /* } */
    maxV = std::max(std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]), maxV);
  }
  t += timeStepSize;
  // deallocate all of our memory to avoid leaks
  /* for (int i = 0; i < NumberOfBodies; i++) delete[] forces[i]; */
  /* delete [] forces[0]; */
  delete [] force0;
  delete [] force1;
  delete [] force2;
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
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
