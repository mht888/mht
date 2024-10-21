#include "include.h"
using namespace emp;
using namespace std;

//Initialize random rA, rB
inline void dorArB(PRG prg, BT* r_A, BT* r_B) {
  for (size_t i = 0; i < BETA; i++) {
    prg.random_data(&(r_A[i * BUCKETS]), BUCKETS * sizeof(BT));
    prg.random_data(&(r_B[i * BUCKETS]), BUCKETS * sizeof(BT));
  }

  for (size_t i = 0; i < BUCKETS; i++) {
    for (uint64_t j = 0; j < BETA; j++) {
      r_B[i * BETA + j] = r_B[i * BETA + j] % MODULUS;
      r_A[i * BETA + j] = r_A[i * BETA + j] % MODULUS;

      //rBs must have multiplicative inverses
      while (r_B[i * BETA + j] == 0) {
        prg.random_data(&(r_B[i * BETA + j]), sizeof(BT));
        r_B[i * BETA + j] = r_B[i * BETA + j] % MODULUS;
      }

    }
  }
}


void computeTriples(BT* r_A, BT* s_A, BT* r_B, BT* s_B) {
  //Use fixed seed for debugging
  PRG prg(fix_key);

  //Initialize sAs
  prg.random_data(s_A, BUCKETS * sizeof(BT));

  //Init rA, rB
  dorArB(prg, r_A, r_B);

  IKNP<NetIO>* ot = new IKNP<NetIO>(ios[0]);
  block* r = NULL;
  bool* b = NULL;
  setup_ot<OT<NetIO>>(ot, ios[0], party);
  if (party == ALICE)
  {
      for (size_t i = 0; i < BUCKETS; i++)
      {
          //Init sA
          if (LOGN > 20)
          {
              for (size_t i = 0; i < (1 << (LOGN - 20)); i++) {
                  cout << "i: " << i << endl;
                  computeShares<OT<NetIO>>(ot, ios[0], &(sA[i * (1 << 20)]), &(rA[i * (1 << 20) * BETA]), party, NULL, NULL, TWO_TO_TWENTY_BUCKETS);
              }
          }
          else {
              computeShares<OT<NetIO>>(ot, ios[0], sA, rA, party, NULL, NULL, BUCKETS);
          }
      }
  }
      else
      {
          //Compute sB
          #pragma omp parallel for
          if (LOGN > 20) {
              for (size_t i = 0; i < (1 << (LOGN - 20)); i++) {
                  computeShares<OT<NetIO>>(ot, ios[0], NULL, NULL, party, &(rB[i * (1 << 20) * BETA]), &(sB[i * (1 << 20) * BETA]), TWO_TO_TWENTY_BUCKETS);
              }
          }
          else {
              computeShares<OT<NetIO>>(ot, ios[0], NULL, NULL, party, rB, sB, BUCKETS);
          }
      }
  if (party == ALICE) {
      free(sA);
  }
  else {//BOB
      delete[] r;
      delete[] b;
      free(sB);
  }
  delete ot;
}

inline void bobsOperations(BT* bobTable, BT* fromAlice, BT* r_B, BT* s_B, BT* INV) {
  //We compute result = ((c + h) + sB) / rB

  #pragma omp parallel for
  for (size_t i = 0; i < BUCKETS; i++) {
    for (size_t j = 0; j < BETA; j++) {
      bobTable[i * BETA + j] = ((fromAlice[i] + bobTable[i * BETA + j] + s_B[i * BETA + j]) * INV[r_B[i * BETA + j]] ) % MODULUS;
    }
  }

}
