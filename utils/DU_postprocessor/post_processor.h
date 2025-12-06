#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

class PostProcessor {
private:
  std::ifstream ifs;

public:
  PostProcessor(const std::string &file_name);
  void graphics_u(double h, double T, const std::string &name, int plot_steps);
  void graphics_elliptic(double h, const std::string &name);
  std::vector<double> error(const std::vector<std::vector<double>> &calc,
                            const std::vector<std::vector<double>> &exact);
};

inline PostProcessor::PostProcessor(const std::string &filename)
    : ifs(filename) {
  if (!ifs.is_open()) {
    throw std::runtime_error("can't open file");
  }
}

inline void PostProcessor::graphics_u(double h, double tau,
                                      const std::string &name, int plot_steps) {
  int K, N;
  ifs >> K >> N; // <-- проверь порядок с файлом!

  std::vector<std::vector<double>> u(K, std::vector<double>(N, 0.0));
  for (int i = 0; i < K; ++i)
    for (int j = 0; j < N; ++j)
      ifs >> u[i][j];

  std::ifstream analytic_ifs("exact.txt");
  analytic_ifs >> K >> N;
  std::vector<std::vector<double>> u_exact(K, std::vector<double>(N, 0.0));
  for (int i = 0; i < K; ++i)
    for (int j = 0; j < N; ++j)
      analytic_ifs >> u_exact[i][j];

  std::vector<double> err = error(u, u_exact);
  double max_err_value = *std::max_element(err.begin(), err.end());
  std::cout << std::format("Max error for method {}: {:.6f}\n", name,
                           max_err_value);

  std::vector<int> time_indices = {
      std::max(0, int(K * 0.05)), std::max(0, int(K * 0.25)),
      std::max(0, int(K * 0.50)), std::max(0, int(K * 0.75)),
      std::max(0, int(K * 0.95))};

  std::ofstream temp("temp_data.txt");
  for (int j = 0; j < u[0].size(); ++j) {
    temp << j * h;
    for (int idx : time_indices) {
      temp << " " << u[idx][j];
    }
    temp << "\n";
  }
  temp.close();

  std::ofstream temp_exact("temp_exact.txt");
  for (int j = 0; j < u_exact[0].size(); ++j) {
    temp_exact << j * h;
    for (int idx : time_indices) {
      temp_exact << " " << u_exact[idx][j];
    }
    temp_exact << "\n";
  }
  temp_exact.close();

  std::ofstream temp_error("temp_error.txt");
  for (int i = 0; i < err.size(); ++i) {
    temp_error << i * tau << ' ' << err[i] << '\n';
  }
  std::cout << std::setprecision(5) << std::fixed << std::endl;
  std::cout << err[0] << ' ' << err[1] << std::endl;
  temp_error.close();

  std::stringstream gnuplotCmd;
  gnuplotCmd
      << "set terminal pngcairo enhanced size 1600,800\n"
      << "set output 'solution_and_error.png'\n"
      << "set multiplot layout 1,2\n"

      // Первый график (решение)
      << "set title '" << name << " over time'\n"
      << "set xlabel 'x'\n"
      << "set ylabel 'u(x,t)'\n"
      << "set grid\n"
      << "plot 'temp_data.txt' using 1:2 title 't=0.05T' with lines dashtype "
         "2,\\\n"
      << "     'temp_data.txt' using 1:3 title 't=0.25T' with lines dashtype "
         "2,\\\n"
      << "     'temp_data.txt' using 1:4 title 't=0.50T' with lines dashtype "
         "2,\\\n"
      << "     'temp_data.txt' using 1:5 title 't=0.75T' with lines dashtype "
         "2,\\\n"
      << "     'temp_data.txt' using 1:6 title 't=0.95T' with lines dashtype 2";

  gnuplotCmd
      << ",\\\n     'temp_exact.txt' using 1:2 title 'exact t=0.05T' with "
         "lines,\\\n"
      << "     'temp_exact.txt' using 1:3 title 'exact t=0.25T' with lines,\\\n"
      << "     'temp_exact.txt' using 1:4 title 'exact t=0.50T' with lines,\\\n"
      << "     'temp_exact.txt' using 1:5 title 'exact t=0.75T' with lines,\\\n"
      << "     'temp_exact.txt' using 1:6 title 'exact t=0.95T' with lines";

  // Второй график (ошибки)
  gnuplotCmd
      << "\n\nset title 'Error over time'\n"
      << "set xlabel 't'\n"
      << "set ylabel 'max error'\n"
      << "set grid\n"
      << "plot 'temp_error.txt' using 1:2 with linespoints pt 7 title 'error'\n"
      << "unset multiplot\n";

  std::ofstream script("gnuplot_script.gp");
  script << gnuplotCmd.str();
  script.close();

  system("gnuplot gnuplot_script.gp");
}

// Замените graphics_u на этот метод
inline void PostProcessor::graphics_elliptic(double h,
                                             const std::string &name) {
  int Nx, Ny; // Вместо K, N используем Nx, Ny (размеры сетки)

  // 1. Читаем численное решение
  ifs >> Nx >> Ny;
  std::vector<std::vector<double>> u(Nx, std::vector<double>(Ny));
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      ifs >> u[i][j];

  // 2. Читаем точное решение
  std::ifstream analytic_ifs("exact.txt");
  // Проверка на открытие файла
  if (!analytic_ifs.is_open()) {
    std::cerr << "Error: exact.txt not found!" << std::endl;
    return;
  }

  int Nx_ex, Ny_ex;
  analytic_ifs >> Nx_ex >> Ny_ex;
  std::vector<std::vector<double>> u_exact(Nx_ex, std::vector<double>(Ny_ex));
  for (int i = 0; i < Nx_ex; ++i)
    for (int j = 0; j < Ny_ex; ++j)
      analytic_ifs >> u_exact[i][j];

  std::ofstream data_file("plot_data.txt");
  double max_diff = 0.0;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      double x = (i + 1) * h;
      double y = (j + 1) * h;
      double diff = std::abs(u[i][j] - u_exact[i][j]);
      if (diff > max_diff)
        max_diff = diff;

      // x y u_num u_exact diff
      data_file << x << " " << y << " " << u[i][j] << " " << u_exact[i][j]
                << " " << diff << "\n";
    }
    data_file << "\n"; // Пустая строка для изоляции рядов (нужно для splot)
  }
  data_file.close();

  std::cout << "Max error between exact solution (" << name << "): " << max_diff
            << std::endl;

  // 4. Скрипт Gnuplot
  std::stringstream gp;
  gp << "set terminal pngcairo enhanced size 1200,600\n"
     << "set output 'elliptic_result.png'\n"
     << "set multiplot layout 1,2 title 'Elliptic Solution: " << name << "'\n"
     << "set title 'Solution Surfaces'\n"
     << "set xlabel 'x'\n"
     << "set ylabel 'y'\n"
     << "set zlabel 'u(x,y)'\n"
     << "set grid\n"
     << "set view 60, 30\n" // Угол обзора
     << "splot 'plot_data.txt' using 1:2:4 with lines title 'Exact' lc rgb "
        "'black', \\\n"
     << "      'plot_data.txt' using 1:2:3 with pm3d title 'Numerical'\n"

     << "set title 'Error Heatmap'\n"
     << "set view map\n" // Вид сверху (2D)
     << "set xlabel 'x'\n"
     << "set ylabel 'y'\n"
     << "set palette defined (0 'white', 1 'red')\n" // Чем краснее, тем больше
                                                     // ошибка
     << "splot 'plot_data.txt' using 1:2:5 with image title 'Error |u-u*|'\n"

     << "unset multiplot\n";

  std::ofstream script("gnuplot_script.gp");
  script << gp.str();
  script.close();

  system("gnuplot gnuplot_script.gp");
}

inline std::vector<double>
PostProcessor::error(const std::vector<std::vector<double>> &calc,
                     const std::vector<std::vector<double>> &exact) {
  size_t K = std::min(calc.size(), exact.size());
  if (K <= 1)
    return {}; // Недостаточно слоев для вычисления ошибки

  std::vector<double> max_errors(K - 1, 0.0);
  for (size_t i = 1; i < K; ++i) { // До i = K-1 (последний слой)
    double max_err = 0.0;
    size_t N = std::min(calc[i].size(), exact[i].size());

    for (size_t j = 0; j < N; ++j) { // От j=0 до N-1
      double err = std::abs(calc[i][j] - exact[i][j]);
      if (err > max_err) {
        max_err = err;
      }
    }
    max_errors[i - 1] = max_err;
  }

  return max_errors;
}