#ifndef PROCESS_BAR_JJ_H
#define PROCESS_BAR_JJ_H

#include <string>
#include <iomanip>
template<typename OS>
void show_progress_bar(OS& os, const double percentage,
                       const std::string &message)
{
  std::string head = message + " ";
  const size_t bar_length = 60;
  if (message.length() >= 10)
    head = std::string(message.begin(), message.begin()+10) + " ";
  
  size_t front_num = std::min(percentage, 100.0) / 100.0 * bar_length;
  size_t back_num = bar_length - front_num;

  os << "\r [" << std::setw(3) << static_cast<int>(percentage) << "%] "
     << head << std::string(front_num, '#') << std::string(back_num, '.') << std::flush;
}

#endif // PROCESS_BAR_JJ_H

