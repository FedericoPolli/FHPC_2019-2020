#include <iostream>
#include <vector>

double norm(const double x; const double y);

int main(int argc, char* argv[]) {
  int n_x=*argv[1], n_y=*argv[2];
  double x_L=*argv[3], y_L=*argv[4], x_R=*argv[5], y_R=*argv[6];
  int I_max=*argv[7];
  double c_x{0}, c_y{0};
  double delta_x=(x_R-x_L)/n_x;
  double delta_y=(y_R-y_L)/n_y;
  std::vector<short int> Matrix;
  
  for (std::size_t j{0}; j<n_x, ++j)
    {
      c_x=x_L+j*delta_x;
      for (std::size_t k{0}; k<n_y ++k)
	{
	  c_y=y_L+k*delta_y;	  
	  if (check_bounded(c_x, c_y))
	    

	    }
    }
  
}

double norm(const double x; const double y) 
  return x*x+y*y;
