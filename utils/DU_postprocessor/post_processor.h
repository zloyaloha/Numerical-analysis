#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>

class PostProcessor {
private:
    std::ifstream ifs;
public:
    PostProcessor(const std::string& file_name);
    void graphics_u(double h, double T, const std::string& name, int plot_steps);
    std::vector<double> error(const std::vector<std::vector<double>>& calc,
                              const std::vector<std::vector<double>>& exact);
};

inline PostProcessor::PostProcessor(const std::string& filename) : ifs(filename)
{
    std::cout << filename << std::endl;
    if (!ifs.is_open()) {
        throw std::runtime_error("can't open file");
    }
}

inline void PostProcessor::graphics_u(double h,
                                               double tau,
                                               const std::string& name,
                                               int plot_steps)
{
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
    std::cout << std::format("Max error for method {}: {:.6f}\n", name, max_err_value);

    std::vector<int> time_indices = {
        std::max(0, int(K * 0.05)),
        std::max(0, int(K * 0.50)),
        std::max(0, int(K * 0.95))
    };

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
    for (int i = 1; i < err.size(); ++i) {
        temp_error << i * tau << ' ' << err[i] << '\n';
    }
    temp_error.close();

    std::stringstream gnuplotCmd;
    gnuplotCmd << "set terminal pngcairo enhanced size 1600,800\n"
               << "set output 'solution_and_error.png'\n"
               << "set multiplot layout 1,2\n"

               // Первый график (решение)
               << "set title '" << name << " over time'\n"
               << "set xlabel 'x'\n"
               << "set ylabel 'u(x,t)'\n"
               << "set grid\n"
               << "plot 'temp_data.txt' using 1:2 title 't=0.05T' with lines dashtype 2,\\\n"
               << "     'temp_data.txt' using 1:3 title 't=0.10T' with lines dashtype 2,\\\n"
               << "     'temp_data.txt' using 1:4 title 't=0.95T' with lines dashtype 2";

    gnuplotCmd << ",\\\n     'temp_exact.txt' using 1:2 title 'exact t=0.05T' with lines,\\\n"
                << "     'temp_exact.txt' using 1:3 title 'exact t=0.10T' with lines,\\\n"
                << "     'temp_exact.txt' using 1:4 title 'exact t=0.95T' with lines";

    // Второй график (ошибки)
    gnuplotCmd << "\n\nset title 'Error over time'\n"
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

inline std::vector<double> PostProcessor::error(const std::vector<std::vector<double>>& calc,
                                                  const std::vector<std::vector<double>>& exact)
{
    // K - 1: Инициализируем массив для K-1 шагов по времени (от k=1 до K-1).
    // Элемент [0] будет хранить ошибку для k=1.
    size_t K = std::min(calc.size(), exact.size());
    if (K <= 1) return {}; // Недостаточно слоев для вычисления ошибки

    std::vector<double> max_errors(K - 1, 0.0);
    // Внешний цикл: итерация по временным слоям, начиная с i = 1 (t=tau)
    for (size_t i = 1; i < K; ++i) { // До i = K-1 (последний слой)
        double max_err = 0.0;
        size_t N = std::min(calc[i].size(), exact[i].size());

        // Внутренний цикл: итерация по пространственным узлам, включая границы (j=0 и j=N-1)
        for (size_t j = 0; j < N; ++j) { // От j=0 до N-1
            double err = std::abs(calc[i][j] - exact[i][j]);
            if (err > max_err) {
                max_err = err;
            }
        }
        // Запись ошибки слоя i в элемент i-1 массива max_errors
        max_errors[i - 1] = max_err;
    }

    return max_errors;
}