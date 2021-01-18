//
// Created by Роман Агеев on 01/12/20.
//
#include <mutex>
#include <future>
#include "quadratic_sieve.hpp"

namespace {
using vec_u64 = std::vector<uint64_t>;
using vec_sz = std::vector<size_t>;
}

std::optional<mpz_class> check_vector(vec_sz const &number_mask, vec_u64 const &factor_primes,
                                      mpz_class const &sqrt_of_number, mpz_class const &number,
                                      std::vector<vec_sz> const & smooth_vectors, vec_sz const & smooth_numbers) {
  vec_u64 counter(factor_primes.size(), 0);
  mpz_class first_operand = 1;
  mpz_class second_operand = 1;
  for (size_t i = 0; i < smooth_vectors.size(); ++i) {
    if (number_mask[i] == 1) {
      second_operand *= (smooth_numbers[i] + sqrt_of_number);
      for (size_t j : smooth_vectors[i]) {
        ++counter[j];
      }
    }
  }
  for (size_t i = 0; i < factor_primes.size(); ++i) {
    for (size_t j = 0; j < (counter[i] / 2); ++j) {
      first_operand *= (size_t) factor_primes[i];
    }
  }
  if ((first_operand % number != (-second_operand) % number + number)
      && (first_operand % number != second_operand % number)) {
    mpz_class ans;
    mpz_gcd(ans.get_mpz_t(), mpz_class(second_operand - first_operand).get_mpz_t(), number.get_mpz_t());
    return ans;
  }
  return {};
}

vec_sz NumberRepresentation::get_number() const {
//  if (row_idx == 0 && col_idx == 0) {
//    mpz_class one = 1;
//    for (auto &i : now) {
//      one *= smooth_numbers[i];
//    }
//    mpz_class r = one + number;
//    mpz_sqrt(r.get_mpz_t(), r.get_mpz_t());
//    mpz_class x{};
//    mpz_sqrt(x.get_mpz_t(), one.get_mpz_t());
//    mpz_class val1 = x - r;
//    mpz_class val2 = x + r;
//    mpz_gcd(val1.get_mpz_t(), val1.get_mpz_t(), number.get_mpz_t());
//    mpz_gcd(val2.get_mpz_t(), val2.get_mpz_t(), number.get_mpz_t());
//    if (one == 1) {
//      return ans;
//    }
//    if (val1 * val2 == number && val1 != 1 && val2 != 1) {
//      ans = {val1, val2};
//    }
//    return ans;
//  } else if (col_idx == primes.size()) {
  size_t random_number = 1;
  NumberRepresentation tmp_number_matrix(*this);
  vec_sz ans(w() - 1, 0);
  for (size_t row_idx = h(); row_idx--; random_number = rand()) {
    size_t sum = 0;
    size_t last_fixed = 0;
    for (size_t col_idx = 0; col_idx < w() - 1; ++col_idx) {
      last_fixed = (bool) tmp_number_matrix(row_idx, col_idx) ? col_idx : last_fixed;
      sum += (bool) tmp_number_matrix(row_idx, col_idx);
    }
    if (sum != 0) {
      size_t val_in_a_box_row = 0;
      if (sum > 1) {
        val_in_a_box_row = random_number % 2;
      } else {
        val_in_a_box_row = (bool) tmp_number_matrix(row_idx, w() - 1);
      }
      ans[last_fixed] = val_in_a_box_row;
      for (size_t j = 0; j <= row_idx; ++j) {
        if (tmp_number_matrix(j, last_fixed)) {
          if (val_in_a_box_row == 1) {
            tmp_number_matrix(j, w() - 1).not_value();
          }
          tmp_number_matrix(j, last_fixed) = false;
        }
      }
      if (sum != 1) {
        ++row_idx;
      }
    }
  }
  return ans;
}

void NumberRepresentation::get_sparce_nullspace() {
  auto &input = (*this);
  for (size_t now_row = 0, now_col = 0; now_row < h() && now_col < w(); ++now_col) {
    size_t min = now_row;
    for (size_t k = now_row + 1; k < h(); ++k) {
      if (input(k, now_col)) {
        min = k;
        break;
      }
    }
//    std::cout << "now matrix:" << std::endl;
//    for (auto & i : input) {
//      for (auto & j : i) {
//        std::cout << j << " ";
//      }
//      std::cout << std::endl;
//    }
//    std::cout << std::endl;
    if (input(min, now_col)) {
      std::swap(matrix[now_row], matrix[min]);
      for (size_t k = now_row + 1; k < h(); ++k) {
        if (input(k, now_col)) {
          for (size_t i_ = 0; i_ < number_boxes_on_row; ++i_) {
            matrix[k][i_] ^= matrix[now_row][i_];
          }
        }
      }
      ++now_row;
    }
  }
}

vec_u64 prime_array(const mpz_class &number, uint64_t B) {
  vec_u64 ans;
  std::vector<bool> mask(B + 1, false);
  for (size_t i = 2; i <= B; ++i) {
    if (!mask[i]) {
      if (mpz_legendre(number.get_mpz_t(), mpz_class(i).get_mpz_t()) == 1) {
        ans.push_back(i);
      }
      for (size_t j = i; j <= B; j += i) {
        mask[j] = true;
      }
    }
  }
  return ans;
}

uint64_t mod_pow(uint64_t fst, uint64_t exp, uint64_t modulus) {
  uint64_t ans = 1;
  while (exp > 0) {
    if (exp % 2 == 1) {
      ans = (ans * fst) % modulus;
    }
    fst *= fst;
    fst %= modulus;
    exp /= 2;
  }
  return ans;
}

vec_sz num_to_vec(mpz_class & quadratic, vec_u64 const & factor_primes) {
  vec_sz factors;
  for (size_t j = 0; j < factor_primes.size(); ++j) {
    while (mpz_divisible_ui_p(quadratic.get_mpz_t(), factor_primes[j])) {
      mpz_divexact_ui(quadratic.get_mpz_t(), quadratic.get_mpz_t(), factor_primes[j]);
      factors.push_back(j);
    }
  }
  return factors;
}

// algorithm from https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
std::pair<uint64_t, uint64_t> tonelli_shanks(uint64_t number, uint64_t prime) {
  if (prime == 2) {
    return {number, number};
  }
  uint64_t sieving = 0;
  uint64_t z_value = 2;
  uint64_t quadratic = prime - 1;
  while (mod_pow(z_value, (prime - 1) / 2, prime) <= 1) {
    ++z_value;
  }
  while (quadratic % 2 == 0) {
    quadratic /= 2;
    ++sieving;
  }
  uint64_t square_counter = mod_pow(z_value, quadratic, prime);
  uint64_t t_value = mod_pow(number, quadratic, prime);
  uint64_t reduses = mod_pow(number, (quadratic + 1) / 2, prime);
  uint64_t tmp_value = sieving;
  while (t_value % prime != 1) {
    uint64_t i = 1;
    while (mod_pow(t_value, std::pow(2, i), prime) != 1) {
      ++i;
    }
    auto pow = std::pow(2, tmp_value - 1 - i);
    tmp_value = i;
    uint64_t pow_mow_result = mod_pow(square_counter, pow, prime);
    square_counter = pow_mow_result * pow_mow_result % prime;
    reduses = reduses * pow_mow_result % prime;
    t_value = t_value * pow_mow_result * pow_mow_result % prime;
  }
  return {reduses, prime - reduses};
}

uint32_t positive_mod(mpz_class const &val, uint32_t mod) {
  return mpz_class((val % mod) + mod).get_ui();
}

auto indexes(vec_u64 const &factor_primes, mpz_class const &sqrt_of_number, mpz_class const &number) {
  std::pair<vec_u64, vec_u64> ans
      {vec_u64(factor_primes.size()), vec_u64(factor_primes.size())};
  for (size_t i = 0; i < factor_primes.size(); ++i) {
    std::pair<uint64_t, uint64_t>
        result = tonelli_shanks(mpz_class(number % (uint32_t) factor_primes[i]).get_ui(), factor_primes[i]);
    ans.first[i] = positive_mod((uint32_t) result.first - sqrt_of_number, factor_primes[i]);
    ans.second[i] = positive_mod((uint32_t) result.second - sqrt_of_number, factor_primes[i]);
  }
  return ans;
}

auto get_residues(size_t array_size, mpz_class const &sqrt_of_number, mpz_class const &number) {
  double last_bound = 0;
  double prev_bound = 1;
  // to use real numbers
  std::vector<double> residues(array_size, 0);
  size_t j = residues.size() + 1;
  for (size_t i = 1; i < residues.size(); ++i) {
    if (prev_bound <= j) {
      prev_bound = 2 * prev_bound + 1;
      last_bound = mpz_sizeinbase(mpz_class((j + sqrt_of_number) * (j + sqrt_of_number) - number).get_mpz_t(), 2);
    }
    residues[i] = last_bound;
    ++j;
  }
  return residues;
}

uint64_t get_number_size(uint64_t num) {
  uint64_t ans = 0;
  while(num != 0) {
    ++ans;
    num >>= 1;
  }
  return ans;
}

mpz_class QuadraticSieve::factorize(mpz_class const &number) {
  auto factor_primes = prime_array(number, B);
  mpz_class sqrt_of_number = sqrt(number);
  uint64_t start_of_interval = 0;
  uint64_t end_of_interval = sieving_interval;
  vec_sz smooth_numbers;
  auto [tonelli_first, tonelli_second] = indexes(factor_primes, sqrt_of_number, number);
  std::vector<vec_sz> smooth_vectors;
//  while (x <= num) {
//    mpz_class opp = x * x - num;
//    if (is_smooth(opp, primes)) {
//      smooth_vectors.push_back(std::move(num_to_vec(opp, primes)));
//      smooth_numbers.push_back(std::move(opp));
//    }
//    x += 1;
//    if (smooth_vectors.size() == n) {
//      smooth_vectors = get_sparce_nullspace(transpose(std::move(smooth_vectors), primes.size()));
//      result = get_number(smooth_numbers, smooth_vectors, primes, num);
//      smooth_vectors.clear();
//      smooth_numbers.clear();
//      if (result) {
//        return result;
//      }
//    }
//  }
  while (factor_primes.size() + min_number > smooth_numbers.size()) {
    auto residues_array = get_residues(sieving_interval, sqrt_of_number, number);
    for (size_t i = 0; i < factor_primes.size(); ++i) {
      auto while_tonelli = [&residues_array, i, &factor_primes, start_of_interval, end_of_interval](vec_u64 &tonelli) {
        while (tonelli[i] < end_of_interval) {
          if (tonelli[i] >= start_of_interval) {
            residues_array[tonelli[i] - start_of_interval] -= get_number_size(factor_primes[i]);
            tonelli[i] += factor_primes[i];
          }
        }
      };
      while_tonelli(tonelli_first);
      while_tonelli(tonelli_second);
    }
    auto last_size = get_number_size(factor_primes.back());
    // can be done in multiple threads
    /*
    auto mutex = std::mutex();
    std::vector<std::future<void>> tasks;
    auto const & residues_array_const = residues_array;
    auto const & sqrt_of_number_const = sqrt_of_number;
    auto const & number_const = number;
    auto const & factor_primes_const = factor_primes;
    for (size_t i = 0; i < sieving_interval; ++i) {
      tasks.emplace_back(std::move(std::async(std::launch::async,
          [i, last_size, start_of_interval, &residues_array_const, &smooth_numbers,
           &smooth_vectors, &mutex, &sqrt_of_number_const,
              &number_const, &factor_primes_const] () {
        if (last_size > std::abs(residues_array_const[i])) {
          mpz_class quadratic = ((size_t)start_of_interval + i + sqrt_of_number_const) *
              ((size_t)start_of_interval + i + sqrt_of_number_const) - number_const;
          auto num_vector = num_to_vec(quadratic, factor_primes_const);
          if (quadratic == 1) {
            mutex.lock();
            smooth_numbers.push_back(start_of_interval + i);
            smooth_vectors.push_back(num_vector);
            mutex.unlock();
          }
        }
      })));
    }
    */
    // but still need some improves
    for (size_t i = 0; i < sieving_interval; ++i) {
      if (last_size > std::abs(residues_array[i])) {
        mpz_class quadratic = ((size_t)start_of_interval + i + sqrt_of_number) *
            ((size_t)start_of_interval + i + sqrt_of_number) - number;
        auto num_vector = num_to_vec(quadratic, factor_primes);
        if (quadratic == 1) {
          smooth_numbers.push_back(start_of_interval + i);
          smooth_vectors.push_back(num_vector);
        }
      }
    }
    end_of_interval += sieving_interval;
    start_of_interval += sieving_interval;
  }
  auto matrix = NumberRepresentation(factor_primes.size(), smooth_vectors.size() + 1, smooth_vectors);
  matrix.get_sparce_nullspace();
  std::optional<mpz_class> ans = {};
  while (!ans.has_value()) {
    ans = check_vector(matrix.get_number(), factor_primes, sqrt_of_number, number, smooth_vectors, smooth_numbers);
  }
  return ans.value();
}
