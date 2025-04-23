// This file was adapted from https://github.com/pq-crystals/kyber (Public Domain)
// authored by: Joppe Bos, Léo Ducas, Eike Kiltz, Tancrède Lepoint, Vadim Lyubashevsky,
// John Schanck, Peter Schwabe, Gregor Seiler, and Damien Stehlé.

#ifndef CPUCYCLES_H
#define CPUCYCLES_H

#include <stdint.h>

#ifdef USE_RDPMC  /* Needs echo 2 > /sys/devices/cpu/rdpmc */

static inline uint64_t cpucycles(void) {
  const uint32_t ecx = (1U << 30) + 1;
  uint64_t result;

  __asm__ volatile ("rdpmc; shlq $32,%%rdx; orq %%rdx,%%rax"
    : "=a" (result) : "c" (ecx) : "rdx");

  return result;
}

#else

static inline uint64_t cpucycles(void) {
  uint64_t result;

  __asm__ volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
    : "=a" (result) : : "%rdx");

  return result;
}

#endif

uint64_t cpucycles_overhead(void);

#endif

void print_results(const char *s, uint64_t *t, size_t tlen);
void get_results(uint64_t *average_value, uint64_t *median_value, uint64_t *t, size_t tlen);
