#include <iostream>
#include <vector>

double check_bounded(const double c_x, const double c_y, const int I_max);
void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
void write_ppm_image( void *image, int maxval, int xsize, int ysize, const char *image_name);
 

int main(int argc, char* argv[]) {
  
  int n_x=0, n_y=0;
  double x_L=0.0, y_L=0.0, x_R=0.0, y_R=0.0;
  int I_max=0;
  
  if (argc>7) {
    n_x=std::stoi(argv[1]);
    n_y=std::stoi(argv[2]);
    x_L=std::stod(argv[3]);
    y_L=std::stod(argv[4]);
    x_R=std::stod(argv[5]);
    y_R=std::stod(argv[6]);
    I_max=std::stoi(argv[7]);
  }
  else {
    std::cout << "Not enough arguments were given.\n";
     return 0;
  }
  double c_x{0}, c_y{0};
  double delta_x=(x_R-x_L)/(double)n_x;
  double delta_y=(y_R-y_L)/(double)n_y;
  std::vector<char> Matrix_255;
  std::vector<short int> Matrix_65535;
  
  if (I_max>=65536) {
    std::cout << "The number of iterations is too big.\n" << "It should be less than 65536.\n";
    return 0;
  }

  if (I_max<256) {
     Matrix_255.resize(n_x*n_y);
#pragma omp parallel for private(c_x, c_y)
    for (std::size_t j=0; j<n_y; ++j)
      {
	c_y=y_L+j*delta_y;
	for (std::size_t k=0; k<n_x; ++k)
	  {
	    c_x=x_L+k*delta_x;	  
	    Matrix_255[k+j*n_x]=check_bounded(c_x, c_y, I_max);
	  }
      }
     write_pgm_image(Matrix_255.data(), I_max, n_x, n_y, "image.pgm" );
  }
  else {
    Matrix_65535.resize(n_x*n_y);
    #pragma omp parallel for private(c_x, c_y)
    for (std::size_t j=0; j<n_y; ++j)
      {
	c_y=y_L+j*delta_y;
	for (std::size_t k=0; k<n_x; ++k)
	  {
	    c_x=x_L+k*delta_x;	  
	    Matrix_65535[k+j*n_x]=check_bounded(c_x, c_y, I_max);
	  }
      }
     write_ppm_image(Matrix_65535.data(), I_max, n_x, n_y, "image.ppm" );
  }
}







double check_bounded(const double c_x, const double c_y, const int I_max) {
  double x_temp=c_x;
  double y_temp=c_y;
  double mul_x=0;
  double mul_y=0;
  for (std::size_t j=1; j<=I_max; ++j) {
    if ((x_temp*x_temp+y_temp*y_temp)>2)  //norm
      return j;
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
     PPM   P3     P6       .ppm        [0-2^16[
  
     ------------------------------------------------------------------ */
}


void write_ppm_image( void *image, int maxval, int xsize, int ysize, const char *image_name)
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

  fprintf(image_file, "P6\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite( image, color_depth, xsize*ysize, image_file);  

  fclose(image_file); 
  return ;

  /* ---------------------------------------------------------------

     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
     ASCII  BINARY

     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16[
  
     ------------------------------------------------------------------ */
}



