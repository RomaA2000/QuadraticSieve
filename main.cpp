#include <iostream>
#include <gmpxx.h>
#include <cmath>
#include <optional>
#include <chrono>
#include "quadratic_sieve.hpp"

uint64_t get_B(mpz_class const & number) {
  const double log_n = mpz_sizeinbase(number.get_mpz_t(), 2);
  const double log_log_n = std::log(log_n);
  return 300 + std::ceil(std::exp(0.5 * std::sqrt(log_n * log_log_n)));
}

void test(mpz_class const & num) {
  auto B = get_B(num);
  auto qs = QuadraticSieve(B, 50, 100000);
  mpz_class d = qs.factorize(num);
  mpz_class r, q;
  mpz_cdiv_q(q.get_mpz_t(), num.get_mpz_t(), d.get_mpz_t());
  mpz_cdiv_r(r.get_mpz_t(), num.get_mpz_t(), d.get_mpz_t());
  mpz_class result = d * q;
  std::cout << "test " << result.get_str() << " = " << d.get_str() << " * " << q.get_str() << std::endl;
  if ((r != 0) || (result != num)) {
    std::cout << "test failed" << std::endl;
  }
}

void test_time(mpz_class const & num) {
  auto B = get_B(num);
  auto qs = QuadraticSieve(B, 50, 100000);
  auto start = std::chrono::steady_clock::now();
  auto repeat_number = 300;
  double sum = 0;
  mpz_class d;
  for (size_t i = 0; i < repeat_number; ++i) {
    d = qs.factorize(num);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    sum += elapsed_seconds.count();
  }
  std::cout << num.get_str() << " " << sum / repeat_number << " " << B << std::endl;
  mpz_class r, q;
  mpz_cdiv_q(q.get_mpz_t(), num.get_mpz_t(), d.get_mpz_t());
  mpz_cdiv_r(r.get_mpz_t(), num.get_mpz_t(), d.get_mpz_t());
  mpz_class result = d * q;
  //std::cout << "test " << result.get_str() << " = " << d.get_str() << " * " << q.get_str() << std::endl;
  assert(r == 0);
  assert(result == num);
}

std::string repeat(size_t size, size_t num) {
  std::string ans = "";
  for (size_t i = 0; i < size; ++i) {
    ans += std::to_string(num);
  }
  return ans;
}

void run_tests() {
  test(1000);
  test(10001);
  test(999999);
  test(2138212);
  test(9999999999);
  test(13472847212);
  test(41573857473);
  test(13472847212);
  test(27847382740);
  test(1232321908);
  test(328348273498);
  test(mpz_class(999979) * mpz_class(999983));
  test(321489240192481);
  test(998381923828938102);
  test(1238917239711238791);
  test(mpz_class("99999999999999999999999999"));
  test(mpz_class("9999999999999999999999999999"));
  test(mpz_class("9182317283712837192837824621"));
  test(mpz_class("987658695847362574856378373839"));
}

int main() {

//  run_tests();

  std::cout << "Введите число:" << std::endl;
  std::string input;
  std::cin >> input;
  test(mpz_class(input));

//  freopen("results.txt", "w", stdout);
//  test_time(2138212);
//  test_time(5232673);
//  test_time(23127382);
//  test_time(273827332);
//  test_time(231273823);
//  test_time(328372872);
//  test_time(2367267326);
//  test_time(9813292839);
//  test_time(21378927154);
//  test_time(98231748135);
//  test_time(271384728132);
//  test_time(982379482798);
//  test_time(1823479827134);
//  test_time(3981294839126);
//  test_time(21340982139048);
//  test_time(99821748231734);
  return 0;
}