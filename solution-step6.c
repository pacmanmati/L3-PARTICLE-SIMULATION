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
  int numBuckets = 10;
  double vBucket = maxV / numBuckets;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();
  double *force0 = new double[NumberOfBodies]();
  double *force1 = new double[NumberOfBodies]();
  double *force2 = new double[NumberOfBodies]();
  /* double **forces = new double*[NumberOfBodies]; // cleanup forces into one pointer for nicer code */
  // appending empty () to the new operator initialises with 0
  /* for (int i = 0; i < NumberOfBodies; i++) forces[i] = new double[3](); */
  // convention: always arr[i][dim]
  double diameter = 0.2;
  int* particleInBucket = new int[NumberOfBodies]();
  // sort the particles into buckets based on velocity
  for (int b = 1; b < numBuckets+1; b++)
    for (int p = 0; p < NumberOfBodies; p++) {
      // if we've assigned a bucket, skip this particle
      if (particleInBucket[p]) continue;
      else if (v[p][0]*v[p][0] + v[p][1]*v[p][1] + v[p][2]*v[p][2] <
	       (b+1)*vBucket*vBucket) particleInBucket[p] = b;
    }
  // place all unsorted particles into the first bucket
  for (int p = 0; p < NumberOfBodies; p++) if (!particleInBucket[p]) particleInBucket[p] = 1;
  for (int bucket = 1; bucket < numBuckets+1; bucket++) {
    int j = 0;
    double oldTimeStepSize = timeStepSize;
    timeStepSize /= std::pow(2, bucket - 1);
    // split dt into particleInBucket[j] portions, then perform them
    for (int iter = 0; iter < std::pow(2, bucket-1); iter++) {
      // find all the particles belonging to the bucket
#pragma omp parallel for reduction(+:force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies]) reduction(min:minDx) schedule(dynamic, 1)
      for (int particle = j; particle < NumberOfBodies; particle++) {
	if (particleInBucket[particle] != bucket) {
	  continue;
	}
	j = particle;
	// zero force every partial timestep
	/* for (int dim = 0; dim < 3; dim++) forces[j][dim] = 0; */
	force0[j] = 0;
	force1[j] = 0;
	force2[j] = 0;
	for (int i = 0; i < NumberOfBodies; i++) {
	  if (i == j) continue;
	  double distanceSquared = 0;
	  // find distance between the jth and ith particle
	  for (int dim = 0; dim < 3; dim++) {
	    distanceSquared += (x[j][dim]-x[i][dim]) * (x[j][dim]-x[i][dim]);
	  }
	  double distance = std::sqrt(distanceSquared);
	  // x,y,z force exerted on j by i
	  /* for (int dim = 0; dim < 3; dim++) */
	  /*   forces[j][dim] += (x[i][dim]-x[j][dim]) * mass[i]*mass[j] / (distanceSquared * distance); */
	  force0[j] -= (x[j][0]-x[i][0]) * mass[i]*mass[j] / distance / distance / distance;
	  force1[j] -= (x[j][1]-x[i][1]) * mass[i]*mass[j] / distance / distance / distance;
	  force2[j] -= (x[j][2]-x[i][2]) * mass[i]*mass[j] / distance / distance / distance;
	  // save minDx
	  minDx = std::min(minDx, distance);
	}
	for (int dim = 0; dim < 3; dim++) {
	  x[j][dim] += timeStepSize * v[j][dim];	  
	}
	v[j][0] += timeStepSize * force0[j] / mass[j];
	v[j][1] += timeStepSize * force1[j] / mass[j];
	v[j][2] += timeStepSize * force2[j] / mass[j];
	maxV = std::max(std::sqrt(v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2]), maxV);
	
      }
      for (int particle = j; particle < NumberOfBodies; particle++) {
	if (particleInBucket[particle] != bucket) {
	  continue;
	}
	for (int i = 0; i < NumberOfBodies; i++) {
	  if (i == j) continue;
	  j = particle;
	  double distanceSquared = 0;
	  for (int dim = 0; dim < 3; dim++) {
	    distanceSquared += (x[j][dim]-x[i][dim]) * (x[j][dim]-x[i][dim]);
	  }
	  // check for collisions
	  if (distanceSquared < diameter*diameter) {
	    // conserve momentum
	    for (int dim = 0; dim < 3; dim++) {
	      v[j][dim] = v[i][dim] * mass[i] / (mass[i]+mass[j]) + v[j][dim] * mass[j] / (mass[i]+mass[j]); // change in velocity as a result of fusion
	      x[j][dim] = x[i][dim] * mass[i] / (mass[i]+mass[j]) + x[j][dim] * mass[j] / (mass[i]+mass[j]); // update the position (somewhere in the middle, depends on masses)
	    }
	    mass[i] += mass[j]; // sum masses of collidants (is that a word? it is now)
	    // overwrite last collidant with the (NumberOfBodies-1)th particle and trim redundant particle
	    const int end = --NumberOfBodies;
	    // if we're on the last particle
	    if (NumberOfBodies < 2) {
	      std::cout << "last body @ (X, Y, Z) : (" << x[0][0] << ", " << x[0][1] << ", "<< x[0][2] << ")" << std::endl;
	      /* tFinal = 0; */
	      std::exit(0);
	    }
	    for (int dim = 0; dim < 3; dim++) {
	      x[i][dim] = x[end][dim];
	      v[i][dim] = v[end][dim];
	    }
	    force0[i] = force0[end];
	    force1[i] = force1[end];
	    force2[i] = force2[end];
	    mass[i] = mass[end];
	    // go back an iter to avoid overlooking a particle
	    i--;
	  }
	}
      }
    }
    timeStepSize = oldTimeStepSize; // revert to old timestepsize
  }
  

  t += timeStepSize;
  // deallocate all of our memory to avoid leaks
  /* for (int dim = 0; dim < 3; dim++) delete[] forces[dim]; */
  /* delete[] forces; */
  delete [] force0;
  delete [] force1;
  delete [] force2;
  delete[] particleInBucket;
}


/**
 * main routine.
 *
 * not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
	      << "  snapshot        interval after how many time units to plot. use 0 to switch off plotting" << std::endl
	      << "  final-time      simulated time (greater 0)" << std::endl
	      << "  dt              time step size (greater 0)" << std::endl
	      << std::endl
	      << "examples:" << std::endl
	      << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t one body moving form the coordinate system's centre along x axis with speed 1" << std::endl
	      << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t one spiralling around the other one" << std::endl
	      << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t three body setup from first lecture" << std::endl
	      << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t five body setup" << std::endl
	      << std::endl
	      << "in this naive code, only the first body moves" << std::endl;

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
