//
// Created by Роман Агеев on 01/12/20.
//

#ifndef QUADRATIC_SIEVE_QUADRATIC_SIEVE_HPP
#define QUADRATIC_SIEVE_QUADRATIC_SIEVE_HPP

#include <cstdint>
#include <vector>
#include <gmpxx.h>
#include <cmath>
#include <optional>

class NumberRepresentation {
  typedef uint64_t box;
  static const uint64_t digits_number = std::numeric_limits<box>::digits;

 public:

  class Helper {
   public:
    friend class NumberRepresentation;
    Helper(NumberRepresentation &matrix, uint32_t h, uint32_t w) :
        box_with_value(matrix.matrix[h][w / digits_number]),
        box_position_mask(box(1) << (w % digits_number)) {}
    Helper(Helper const &) = default;
    Helper(Helper &&) = default;
    Helper &operator=(Helper const &) = delete;
    Helper &operator=(Helper &&) = delete;
    Helper &operator=(bool value) {
      if (value) {
        box_with_value = box_with_value | box_position_mask;
      } else {
        box_with_value = box_with_value & ~box_position_mask;
      }
      return *this;
    }

    void not_value() {
      box_with_value ^= box_position_mask;
    }

    explicit operator bool() const {
      return (box_with_value & box_position_mask) != 0;
    }

    box &box_with_value;
    box box_position_mask;
  };

  NumberRepresentation(size_t h, size_t w) : number_boxes_on_row(w / digits_number + 1),
                                             matrix(h, std::vector<box>(w)) {}

  NumberRepresentation(size_t h, size_t w, std::vector<std::vector<size_t>> const &values) : NumberRepresentation(h, w) {
    for (size_t i = 0; i < values.size(); ++i) {
      for (size_t j = 0; j < values[i].size(); ++j) {
        (*this)(values[i][j], i).not_value();
      }
    }
  }

  NumberRepresentation(NumberRepresentation const &) = default;
  NumberRepresentation(NumberRepresentation &&) = default;
  NumberRepresentation & operator=(NumberRepresentation const &) = default;
  NumberRepresentation & operator=(NumberRepresentation &&) = default;
  ~NumberRepresentation() = default;

  Helper operator()(uint32_t h, uint32_t w) const {
    return Helper(const_cast<NumberRepresentation &>(*this), h, w);
  }

  Helper operator()(uint32_t h, uint32_t w) {
    return Helper(*this, h, w);
  }

  void get_sparce_nullspace();

  [[nodiscard]] std::vector<size_t> get_number() const;

  [[nodiscard]] size_t h() const { return matrix.size(); }
  [[nodiscard]] size_t w() const { return matrix[0].size();; }
 private:
  std::vector<std::vector<box>> matrix{};
  size_t number_boxes_on_row = 0;
};

class QuadraticSieve {
 private:
  uint64_t const B;
  uint64_t const min_number;
  uint64_t const sieving_interval;
 public:
  QuadraticSieve(uint64_t B, uint64_t min_number, uint64_t sieving_interval)
      : B(B), min_number(min_number), sieving_interval(sieving_interval) {}

  mpz_class factorize(const mpz_class& number);
};

#endif //QUADRATIC_SIEVE_QUADRATIC_SIEVE_HPP
