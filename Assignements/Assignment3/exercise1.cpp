#include <iostream>
#include <vector>

double norm(const double x, const double y);
auto complex_mul(const double x_1, const double y_1, const double x_2, const double y_2);
double check_bounded(const double c_x, const double c_y, const int I_max);

int main(int argc, char* argv[]) {
  int n_x=std::stoi(argv[1]), n_y=std::stoi(argv[2]);
  double x_L=std::stod(argv[3]), y_L=std::stod(argv[4]), x_R=std::stod(argv[5]), y_R=std::stod(argv[6]);
  int I_max=std::stoi(argv[7]);
  double c_x{0}, c_y{0};
  double delta_x=(x_R-x_L)/n_x;
  double delta_y=(y_R-y_L)/n_y;
  std::vector<short int> Matrix;
  Matrix.reserve(n_x*n_y);
  
  for (std::size_t j{0}; j<n_x; ++j)
    {
      c_x=x_L+j*delta_x;
      for (std::size_t k{0}; k<n_y; ++k)
	{
	  c_y=y_L+k*delta_y;	  
	  Matrix[k+j*n_y]=check_bounded(c_x, c_y, I_max);
	  std::cout << Matrix[k+j*n_y] << " ";
	}
      std::cout << "\n\n";
    }
  
}




double norm(const double x, const double y) {
  return x*x+y*y;
}

auto complex_mul(const double x_1, const double y_1, const double x_2, const double y_2) {
  return std::vector<double> {x_1*x_2-y_1*y_2, x_1*y_2+x_2*y_1};
}

double check_bounded(const double c_x, const double c_y, const int I_max) {
  double x_temp=c_x;
  double y_temp=c_y;
  for (std::size_t j{1}; j<=I_max; ++j) {
    if (norm(x_temp, y_temp)>2)
      return j;
    auto mul=complex_mul(x_temp, y_temp, x_temp, y_temp);
    x_temp=mul[0]+c_x;
    y_temp=mul[1]+c_y;
  }
  return 0;
}
