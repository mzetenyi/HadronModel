#ifndef TABULATOR_HPP
#define TABULATOR_HPP

#include <iomanip>
#include <iostream>
#include <vector>

/**
 * @brief Print tabulated output with predefined field widths.
 * 
 * Example usage: 
 * 
 *     Tabulator tab(5,12);
 *     tab.printComment("sqrt(s)","sigtot");
 *     tab.printComment("[GeV]","[mb]");
 *     for (double srt(srtmin); srt < srtmax; srt += dsrt) {
 *        tab.printLine(GeV(srt),mb(sigtot));
 *     } 
 * 
 */
class Tabulator {
 public:
  template <typename... Args>
  Tabulator(Args... args) {
    (sizes.push_back(args), ...);
  }

  template <typename... Args>
  void printLine(std::ostream& out, Args... args) {
    out << " ";
    printLine(0, out, args...);
  }

  template <typename... Args>
  void printLine(Args... args) {
    cout << " ";
    printLine(0, cout, args...);
  }

  template <typename... Args>
  void printComment(std::ostream& out, Args... args) {
    out << "#";
    printLine(0, out, args...);
  }

  template <typename... Args>
  void printComment(Args... args) {
    cout << "#";
    printLine(0, cout, args...);
  }

private:
  template <typename T, typename... Args>
  void printLine(int i, std::ostream& out, T a1, Args... args) {
    out << " " << std::setw(sizes[i]) << a1;
    printLine(i + 1, out, args...);
  }

  template <typename T>
  void printLine(int i, std::ostream& out, T a1) {
    out << " " << std::setw(sizes[i]) << a1 << endl;
  }

  std::vector<int> sizes;
};

#endif  // TABULATOR_HPP
