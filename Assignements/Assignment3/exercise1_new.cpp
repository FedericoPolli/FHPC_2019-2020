#include <iostream>
#include <vector>
#include <time.h>
#include <omp.h>
#include <exception>
#include <cmath>

#if defined(_OPENMP)
#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec + \
		     (double)myts.tv_nsec * 1e-9)

#endif

double check_bounded(const double c_x, const double c_y, const int I_max);
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);



int main(int argc, char* argv[]) {
  
  int n_x=0, n_y=0;
  double x_L=0.0, y_L=0.0, x_R=0.0, y_R=0.0;
  int I_max=0;

  // Verify that the correct amount of command line arguments were given
  // and that all variables have their appropriate type.
  
  if (argc>7) {
    try {
      n_x=std::stoi(argv[1]);
      n_y=std::stoi(argv[2]);
      x_L=std::stod(argv[3]);
      y_L=std::stod(argv[4]);
      x_R=std::stod(argv[5]);
      y_R=std::stod(argv[6]);
      I_max=std::stoi(argv[7]);
    }
    catch(...) {
      std::cout << "An exception occured, there was probably a typing error.\n\n";
      std::cout << "Terminating...\n\n";
      return 0;
    }
  }
  else {
    std::cout << "Not enough arguments were given.\n\n";
    std::cout << "Input should be n_x, n_y, x_L, y_L, x_R, y_R, I_max, where:\n\n";
    std::cout << "n_x, n_y are the dimensions of the grid of pixels \n\n";
    std::cout << "x_L, y_L, x_R, y_R are the coordinates of the bottom left and top right points \n\n";
    std::cout << "I_max is the number of iterations (less than 65536).\n\n";
    return 0;
  }

  if (I_max != 65535 && I_max !=255 ) {
    std::cout << "Wrong number of iterations: it should be either 255 or 65535.\n\n";
    return 0;
  }

  double c_x{0}, c_y{0};
  double delta_x=(x_R-x_L)/(double)n_x;
  double delta_y=(y_R-y_L)/(double)n_y;
  
  std::vector<char> Matrix_255;
  std::vector<short int> Matrix_65535;
  
  int nthreads=1;
  double avg_time=0;
  double min_time=1e11;
  
  // Distinguishes two cases depending on the number of iterations, which determines what vector to use.
  
  if (I_max==255)
    {
      Matrix_255.resize(n_x*n_y);

#pragma omp parallel private(c_x, c_y) reduction (+:avg_time) reduction (min:min_time) 
      {
	
#if defined(_OPENMP)    
#pragma omp single nowait
	nthreads=omp_get_num_threads();
	
	struct timespec myts;
	double mystart=CPU_TIME_th;
#endif

	
#pragma omp for schedule(dynamic)
	for (std::size_t j=0; j<n_x*n_y; ++j)
	  {
	    c_y=y_L+delta_y*(j/n_x);
	    c_x=x_L+delta_x*(j%n_x);	  
	    Matrix_255[j]=check_bounded(c_x, c_y, I_max);	      
	  }
	
#if defined(_OPENMP)
	avg_time += CPU_TIME_th-mystart;
	min_time = CPU_TIME_th-mystart;
#endif
      }
      write_pgm_image(Matrix_255.data(), I_max, n_x, n_y, "image.pgm" );
    }
  
  else
    {
      Matrix_65535.resize(n_x*n_y);
    
#pragma omp parallel private(c_x, c_y) reduction (+:avg_time) reduction (min:min_time)
      {
      
#if defined(_OPENMP)       
#pragma omp single nowait
	nthreads=omp_get_num_threads();
	
	struct timespec myts;
	double mystart=CPU_TIME_th;
#endif

      
#pragma omp for schedule(dynamic)
	for (std::size_t j=0; j<n_x*n_y; ++j)
	  {
	    c_y=y_L+delta_y*(j/n_x);
	    c_x=x_L+delta_x*(j%n_x);	  
	    Matrix_65535[j]=check_bounded(c_x, c_y, I_max);	      
	  }

#if defined(_OPENMP)
	avg_time += CPU_TIME_th-mystart;
	min_time = CPU_TIME_th-mystart;
#endif      
      }
      write_pgm_image(Matrix_65535.data(), I_max, n_x, n_y, "image.pgm" );
    }

#if defined(_OPENMP)
  std::cout << "Average thread time " << avg_time/nthreads << std::endl;
  std::cout << "Minimum thread time " << min_time << std::endl;
#endif 

}


double check_bounded(const double c_x, const double c_y, const int I_max)
{
  /*------------------------------------------------------------------------------------------
    
    Takes as input the number of iterations and the starting point c=c_1=c_x+i*c_y.

    Taking into account the fact that the Mandelbrot set has two well defined subregions, a 
    cardioid on the right and a disk on the left, the function first checks if the given point
    belongs to eiher the former or the latter, to avoid unnecessary computations when 
    calculating the whole set (gains 2x-3x speedup). Those checks should be commented out if 
    trying to zoom in the frontier.
 
    Then in every iteration it first checks whether |c_n|>2, and if not calculates 
    c_(n+1)=(c_n)^2+c. It either returns 0 if the point belongs to the Mandelbrot set 
    or the iteration at which the norm became greater than 2.    
                             
    ------------------------------------------------------------------------------------------*/

  
  if (sqrt((c_x+0.25)*(c_x+0.25)+c_y*c_y)<0.5)  // Considering the cardioid on the right, the biggest disk
    return 0;                                   // contained in it is centered in (-0.25, 0) with radius 0.5
  
  if (sqrt((c_x+1)*(c_x+1)+c_y*c_y)<0.25)      // The left disk is centered in (-1, 0) with radius 0.25
    return 0;
  
  double x_temp=c_x;
  double y_temp=c_y;
  double mul_x=0;
  double mul_y=0;
  std::size_t j=1;

  for (j=1; j<=I_max; ++j)
    {
      if ((x_temp*x_temp+y_temp*y_temp)>2)  
	return I_max-j;
      mul_x=(x_temp*x_temp)-(y_temp*y_temp);
      mul_y=x_temp*y_temp+x_temp*y_temp;
      x_temp=mul_x+c_x;
      y_temp=mul_y+c_y;
    }
  return 0;   
}


void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name)
{
  FILE* image_file; 
  image_file = fopen(image_name, "w"); 
  
  // Writing header
  // The header's format is as follows, all in ASCII.
  // "whitespace" is either a blank or a TAB or a CF or a LF
  // - The Magic Number (see below the magic numbers)
  // - the image's width
  // - the height
  // - a white space
  // - the image's height
  // - a whitespace
  // - the maximum color value, which must be between 0 and 65535
  //
  // if he maximum color value is in the range [0-255], then
  // a pixel will be expressed by a single byte; if the maximum is
  // larger than 255, then 2 bytes will be needed for each pixel
  //

  int color_depth = 1+((maxval>>8)>0);       // 1 if maxval < 256, 2 otherwise

  fprintf(image_file, "P5\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite( image, color_depth, xsize*ysize, image_file);  

  fclose(image_file); 
  return ;

  /* ---------------------------------------------------------------

     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
     ASCII  BINARY

     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16]
  
     ------------------------------------------------------------------ */
}

