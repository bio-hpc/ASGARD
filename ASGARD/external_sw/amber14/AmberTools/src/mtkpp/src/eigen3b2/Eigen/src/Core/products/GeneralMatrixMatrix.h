// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef EIGEN_GENERAL_MATRIX_MATRIX_H
#define EIGEN_GENERAL_MATRIX_MATRIX_H

template<typename _LhsScalar, typename _RhsScalar> class ei_level3_blocking;

/* Specialization for a row-major destination matrix => simple transposition of the product */
template<
  typename Index,
  typename LhsScalar, int LhsStorageOrder, bool ConjugateLhs,
  typename RhsScalar, int RhsStorageOrder, bool ConjugateRhs>
struct ei_general_matrix_matrix_product<Index,LhsScalar,LhsStorageOrder,ConjugateLhs,RhsScalar,RhsStorageOrder,ConjugateRhs,RowMajor>
{
  typedef typename ei_scalar_product_traits<LhsScalar, RhsScalar>::ReturnType ResScalar;
  static EIGEN_STRONG_INLINE void run(
    Index rows, Index cols, Index depth,
    const LhsScalar* lhs, Index lhsStride,
    const RhsScalar* rhs, Index rhsStride,
    ResScalar* res, Index resStride,
    ResScalar alpha,
    ei_level3_blocking<RhsScalar,LhsScalar>& blocking,
    GemmParallelInfo<Index>* info = 0)
  {
    // transpose the product such that the result is column major
    ei_general_matrix_matrix_product<Index,
      RhsScalar, RhsStorageOrder==RowMajor ? ColMajor : RowMajor, ConjugateRhs,
      LhsScalar, LhsStorageOrder==RowMajor ? ColMajor : RowMajor, ConjugateLhs,
      ColMajor>
    ::run(cols,rows,depth,rhs,rhsStride,lhs,lhsStride,res,resStride,alpha,blocking,info);
  }
};

/*  Specialization for a col-major destination matrix
 *    => Blocking algorithm following Goto's paper */
template<
  typename Index,
  typename LhsScalar, int LhsStorageOrder, bool ConjugateLhs,
  typename RhsScalar, int RhsStorageOrder, bool ConjugateRhs>
struct ei_general_matrix_matrix_product<Index,LhsScalar,LhsStorageOrder,ConjugateLhs,RhsScalar,RhsStorageOrder,ConjugateRhs,ColMajor>
{
typedef typename ei_scalar_product_traits<LhsScalar, RhsScalar>::ReturnType ResScalar;
static void run(Index rows, Index cols, Index depth,
  const LhsScalar* _lhs, Index lhsStride,
  const RhsScalar* _rhs, Index rhsStride,
  ResScalar* res, Index resStride,
  ResScalar alpha,
  ei_level3_blocking<LhsScalar,RhsScalar>& blocking,
  GemmParallelInfo<Index>* info = 0)
{
  ei_const_blas_data_mapper<LhsScalar, Index, LhsStorageOrder> lhs(_lhs,lhsStride);
  ei_const_blas_data_mapper<RhsScalar, Index, RhsStorageOrder> rhs(_rhs,rhsStride);

  typedef ei_gebp_traits<LhsScalar,RhsScalar> Traits;

  Index kc = blocking.kc();                 // cache block size along the K direction
  Index mc = std::min(rows,blocking.mc());  // cache block size along the M direction
  //Index nc = blocking.nc(); // cache block size along the N direction

  ei_gemm_pack_lhs<LhsScalar, Index, Traits::mr, Traits::LhsProgress, LhsStorageOrder> pack_lhs;
  ei_gemm_pack_rhs<RhsScalar, Index, Traits::nr, RhsStorageOrder> pack_rhs;
  ei_gebp_kernel<LhsScalar, RhsScalar, Index, Traits::mr, Traits::nr, ConjugateLhs, ConjugateRhs> gebp;

#ifdef EIGEN_HAS_OPENMP
  if(info)
  {
    // this is the parallel version!
    Index tid = omp_get_thread_num();
    Index threads = omp_get_num_threads();

    LhsScalar* blockA = ei_aligned_stack_new(LhsScalar, kc*mc);
    std::size_t sizeW = kc*Traits::WorkSpaceFactor;
    RhsScalar* w = ei_aligned_stack_new(RhsScalar, sizeW);
    RhsScalar* blockB = blocking.blockB();
    ei_internal_assert(blockB!=0);

    // For each horizontal panel of the rhs, and corresponding vertical panel of the lhs...
    for(Index k=0; k<depth; k+=kc)
    {
      const Index actual_kc = std::min(k+kc,depth)-k; // => rows of B', and cols of the A'

      // In order to reduce the chance that a thread has to wait for the other,
      // let's start by packing A'.
      pack_lhs(blockA, &lhs(0,k), lhsStride, actual_kc, mc);

      // Pack B_k to B' in a parallel fashion:
      // each thread packs the sub block B_k,j to B'_j where j is the thread id.

      // However, before copying to B'_j, we have to make sure that no other thread is still using it,
      // i.e., we test that info[tid].users equals 0.
      // Then, we set info[tid].users to the number of threads to mark that all other threads are going to use it.
      while(info[tid].users!=0) {}
      info[tid].users += threads;

      pack_rhs(blockB+info[tid].rhs_start*kc, &rhs(k,info[tid].rhs_start), rhsStride, actual_kc, info[tid].rhs_length);

      // Notify the other threads that the part B'_j is ready to go.
      info[tid].sync = k;

      // Computes C_i += A' * B' per B'_j
      for(Index shift=0; shift<threads; ++shift)
      {
        Index j = (tid+shift)%threads;

        // At this point we have to make sure that B'_j has been updated by the thread j,
        // we use testAndSetOrdered to mimic a volatile access.
        // However, no need to wait for the B' part which has been updated by the current thread!
        if(shift>0)
          while(info[j].sync!=k) {}

        gebp(res+info[j].rhs_start*resStride, resStride, blockA, blockB+info[j].rhs_start*kc, mc, actual_kc, info[j].rhs_length, alpha, -1,-1,0,0, w);
      }

      // Then keep going as usual with the remaining A'
      for(Index i=mc; i<rows; i+=mc)
      {
        const Index actual_mc = std::min(i+mc,rows)-i;

        // pack A_i,k to A'
        pack_lhs(blockA, &lhs(i,k), lhsStride, actual_kc, actual_mc);

        // C_i += A' * B'
        gebp(res+i, resStride, blockA, blockB, actual_mc, actual_kc, cols, alpha, -1,-1,0,0, w);
      }

      // Release all the sub blocks B'_j of B' for the current thread,
      // i.e., we simply decrement the number of users by 1
      for(Index j=0; j<threads; ++j)
        #pragma omp atomic
        --(info[j].users);
    }

    ei_aligned_stack_delete(LhsScalar, blockA, kc*mc);
    ei_aligned_stack_delete(RhsScalar, w, sizeW);
  }
  else
#endif // EIGEN_HAS_OPENMP
  {
    EIGEN_UNUSED_VARIABLE(info);

    // this is the sequential version!
    std::size_t sizeA = kc*mc;
    std::size_t sizeB = kc*cols;
    std::size_t sizeW = kc*Traits::WorkSpaceFactor;
    LhsScalar *blockA = blocking.blockA()==0 ? ei_aligned_stack_new(LhsScalar, sizeA) : blocking.blockA();
    RhsScalar *blockB = blocking.blockB()==0 ? ei_aligned_stack_new(RhsScalar, sizeB) : blocking.blockB();
    RhsScalar *blockW = blocking.blockW()==0 ? ei_aligned_stack_new(RhsScalar, sizeW) : blocking.blockW();

    // For each horizontal panel of the rhs, and corresponding panel of the lhs...
    // (==GEMM_VAR1)
    for(Index k2=0; k2<depth; k2+=kc)
    {
      const Index actual_kc = std::min(k2+kc,depth)-k2;

      // OK, here we have selected one horizontal panel of rhs and one vertical panel of lhs.
      // => Pack rhs's panel into a sequential chunk of memory (L2 caching)
      // Note that this panel will be read as many times as the number of blocks in the lhs's
      // vertical panel which is, in practice, a very low number.
      pack_rhs(blockB, &rhs(k2,0), rhsStride, actual_kc, cols);


      // For each mc x kc block of the lhs's vertical panel...
      // (==GEPP_VAR1)
      for(Index i2=0; i2<rows; i2+=mc)
      {
        const Index actual_mc = std::min(i2+mc,rows)-i2;

        // We pack the lhs's block into a sequential chunk of memory (L1 caching)
        // Note that this block will be read a very high number of times, which is equal to the number of
        // micro vertical panel of the large rhs's panel (e.g., cols/4 times).
        pack_lhs(blockA, &lhs(i2,k2), lhsStride, actual_kc, actual_mc);

        // Everything is packed, we can now call the block * panel kernel:
        gebp(res+i2, resStride, blockA, blockB, actual_mc, actual_kc, cols, alpha, -1, -1, 0, 0, blockW);

      }
    }

    if(blocking.blockA()==0) ei_aligned_stack_delete(LhsScalar, blockA, kc*mc);
    if(blocking.blockB()==0) ei_aligned_stack_delete(RhsScalar, blockB, sizeB);
    if(blocking.blockW()==0) ei_aligned_stack_delete(RhsScalar, blockW, sizeW);
  }
}

};

/*********************************************************************************
*  Specialization of GeneralProduct<> for "large" GEMM, i.e.,
*  implementation of the high level wrapper to ei_general_matrix_matrix_product
**********************************************************************************/

template<typename Lhs, typename Rhs>
struct ei_traits<GeneralProduct<Lhs,Rhs,GemmProduct> >
 : ei_traits<ProductBase<GeneralProduct<Lhs,Rhs,GemmProduct>, Lhs, Rhs> >
{};

template<typename Scalar, typename Index, typename Gemm, typename Lhs, typename Rhs, typename Dest, typename BlockingType>
struct ei_gemm_functor
{
  ei_gemm_functor(const Lhs& lhs, const Rhs& rhs, Dest& dest, Scalar actualAlpha,
                  BlockingType& blocking)
    : m_lhs(lhs), m_rhs(rhs), m_dest(dest), m_actualAlpha(actualAlpha), m_blocking(blocking)
  {}

  void initParallelSession() const
  {
    m_blocking.allocateB();
  }

  void operator() (Index row, Index rows, Index col=0, Index cols=-1, GemmParallelInfo<Index>* info=0) const
  {
    if(cols==-1)
      cols = m_rhs.cols();

    Gemm::run(rows, cols, m_lhs.cols(),
              /*(const Scalar*)*/&(m_lhs.const_cast_derived().coeffRef(row,0)), m_lhs.outerStride(),
              /*(const Scalar*)*/&(m_rhs.const_cast_derived().coeffRef(0,col)), m_rhs.outerStride(),
              (Scalar*)&(m_dest.coeffRef(row,col)), m_dest.outerStride(),
              m_actualAlpha, m_blocking, info);
  }

  protected:
    const Lhs& m_lhs;
    const Rhs& m_rhs;
    Dest& m_dest;
    Scalar m_actualAlpha;
    BlockingType& m_blocking;
};

template<int StorageOrder, typename LhsScalar, typename RhsScalar, int MaxRows, int MaxCols, int MaxDepth,
bool FiniteAtCompileTime = MaxRows!=Dynamic && MaxCols!=Dynamic && MaxDepth != Dynamic> class ei_gemm_blocking_space;

template<typename _LhsScalar, typename _RhsScalar>
class ei_level3_blocking
{
    typedef _LhsScalar LhsScalar;
    typedef _RhsScalar RhsScalar;

  protected:
    LhsScalar* m_blockA;
    RhsScalar* m_blockB;
    RhsScalar* m_blockW;

    DenseIndex m_mc;
    DenseIndex m_nc;
    DenseIndex m_kc;

  public:

    ei_level3_blocking()
      : m_blockA(0), m_blockB(0), m_blockW(0), m_mc(0), m_nc(0), m_kc(0)
    {}

    inline DenseIndex mc() const { return m_mc; }
    inline DenseIndex nc() const { return m_nc; }
    inline DenseIndex kc() const { return m_kc; }

    inline LhsScalar* blockA() { return m_blockA; }
    inline RhsScalar* blockB() { return m_blockB; }
    inline RhsScalar* blockW() { return m_blockW; }
};

template<int StorageOrder, typename _LhsScalar, typename _RhsScalar, int MaxRows, int MaxCols, int MaxDepth>
class ei_gemm_blocking_space<StorageOrder,_LhsScalar,_RhsScalar,MaxRows, MaxCols, MaxDepth, true>
  : public ei_level3_blocking<
      typename ei_meta_if<StorageOrder==RowMajor,_RhsScalar,_LhsScalar>::ret,
      typename ei_meta_if<StorageOrder==RowMajor,_LhsScalar,_RhsScalar>::ret>
{
    enum {
      Transpose = StorageOrder==RowMajor,
      ActualRows = Transpose ? MaxCols : MaxRows,
      ActualCols = Transpose ? MaxRows : MaxCols
    };
    typedef typename ei_meta_if<Transpose,_RhsScalar,_LhsScalar>::ret LhsScalar;
    typedef typename ei_meta_if<Transpose,_LhsScalar,_RhsScalar>::ret RhsScalar;
    typedef ei_gebp_traits<LhsScalar,RhsScalar> Traits;
    enum {
      SizeA = ActualRows * MaxDepth,
      SizeB = ActualCols * MaxDepth,
      SizeW = MaxDepth * Traits::WorkSpaceFactor
    };

    EIGEN_ALIGN16 LhsScalar m_staticA[SizeA];
    EIGEN_ALIGN16 RhsScalar m_staticB[SizeB];
    EIGEN_ALIGN16 RhsScalar m_staticW[SizeW];

  public:

    ei_gemm_blocking_space(DenseIndex /*rows*/, DenseIndex /*cols*/, DenseIndex /*depth*/)
    {
      this->m_mc = ActualRows;
      this->m_nc = ActualCols;
      this->m_kc = MaxDepth;
      this->m_blockA = m_staticA;
      this->m_blockB = m_staticB;
      this->m_blockW = m_staticW;
    }

    inline void allocateA() {}
    inline void allocateB() {}
    inline void allocateW() {}
    inline void allocateAll() {}
};

template<int StorageOrder, typename _LhsScalar, typename _RhsScalar, int MaxRows, int MaxCols, int MaxDepth>
class ei_gemm_blocking_space<StorageOrder,_LhsScalar,_RhsScalar,MaxRows, MaxCols, MaxDepth, false>
  : public ei_level3_blocking<
      typename ei_meta_if<StorageOrder==RowMajor,_RhsScalar,_LhsScalar>::ret,
      typename ei_meta_if<StorageOrder==RowMajor,_LhsScalar,_RhsScalar>::ret>
{
    enum {
      Transpose = StorageOrder==RowMajor
    };
    typedef typename ei_meta_if<Transpose,_RhsScalar,_LhsScalar>::ret LhsScalar;
    typedef typename ei_meta_if<Transpose,_LhsScalar,_RhsScalar>::ret RhsScalar;
    typedef ei_gebp_traits<LhsScalar,RhsScalar> Traits;

    DenseIndex m_sizeA;
    DenseIndex m_sizeB;
    DenseIndex m_sizeW;

  public:

    ei_gemm_blocking_space(DenseIndex rows, DenseIndex cols, DenseIndex depth)
    {
      this->m_mc = Transpose ? cols : rows;
      this->m_nc = Transpose ? rows : cols;
      this->m_kc = depth;

      computeProductBlockingSizes<LhsScalar,RhsScalar>(this->m_kc, this->m_mc, this->m_nc);
      m_sizeA = this->m_mc * this->m_kc;
      m_sizeB = this->m_kc * this->m_nc;
      m_sizeW = this->m_kc*Traits::WorkSpaceFactor;
    }

    void allocateA()
    {
      if(this->m_blockA==0)
        this->m_blockA = ei_aligned_new<LhsScalar>(m_sizeA);
    }

    void allocateB()
    {
      if(this->m_blockB==0)
        this->m_blockB = ei_aligned_new<RhsScalar>(m_sizeB);
    }

    void allocateW()
    {
      if(this->m_blockW==0)
        this->m_blockW = ei_aligned_new<RhsScalar>(m_sizeW);
    }

    void allocateAll()
    {
      allocateA();
      allocateB();
      allocateW();
    }

    ~ei_gemm_blocking_space()
    {
      ei_aligned_delete(this->m_blockA, m_sizeA);
      ei_aligned_delete(this->m_blockB, m_sizeB);
      ei_aligned_delete(this->m_blockW, m_sizeW);
    }
};

template<typename Lhs, typename Rhs>
class GeneralProduct<Lhs, Rhs, GemmProduct>
  : public ProductBase<GeneralProduct<Lhs,Rhs,GemmProduct>, Lhs, Rhs>
{
    enum {
      MaxDepthAtCompileTime = EIGEN_SIZE_MIN_PREFER_FIXED(Lhs::MaxColsAtCompileTime,Rhs::MaxRowsAtCompileTime)
    };
  public:
    EIGEN_PRODUCT_PUBLIC_INTERFACE(GeneralProduct)
    
    typedef typename  Lhs::Scalar LhsScalar;
    typedef typename  Rhs::Scalar RhsScalar;
    typedef           Scalar      ResScalar;

    GeneralProduct(const Lhs& lhs, const Rhs& rhs) : Base(lhs,rhs)
    {
      typedef ei_scalar_product_op<LhsScalar,RhsScalar> BinOp;
      EIGEN_CHECK_BINARY_COMPATIBILIY(BinOp,LhsScalar,RhsScalar);
    }

    template<typename Dest> void scaleAndAddTo(Dest& dst, Scalar alpha) const
    {
      ei_assert(dst.rows()==m_lhs.rows() && dst.cols()==m_rhs.cols());

      ActualLhsType lhs = LhsBlasTraits::extract(m_lhs);
      ActualRhsType rhs = RhsBlasTraits::extract(m_rhs);

      Scalar actualAlpha = alpha * LhsBlasTraits::extractScalarFactor(m_lhs)
                                 * RhsBlasTraits::extractScalarFactor(m_rhs);

      typedef ei_gemm_blocking_space<(Dest::Flags&RowMajorBit) ? RowMajor : ColMajor,LhsScalar,RhsScalar,
              Dest::MaxRowsAtCompileTime,Dest::MaxColsAtCompileTime,MaxDepthAtCompileTime> BlockingType;

      typedef ei_gemm_functor<
        Scalar, Index,
        ei_general_matrix_matrix_product<
          Index,
          LhsScalar, (_ActualLhsType::Flags&RowMajorBit) ? RowMajor : ColMajor, bool(LhsBlasTraits::NeedToConjugate),
          RhsScalar, (_ActualRhsType::Flags&RowMajorBit) ? RowMajor : ColMajor, bool(RhsBlasTraits::NeedToConjugate),
          (Dest::Flags&RowMajorBit) ? RowMajor : ColMajor>,
        _ActualLhsType, _ActualRhsType, Dest, BlockingType> GemmFunctor;

      BlockingType blocking(dst.rows(), dst.cols(), lhs.cols());

      ei_parallelize_gemm<(Dest::MaxRowsAtCompileTime>32 || Dest::MaxRowsAtCompileTime==Dynamic)>(GemmFunctor(lhs, rhs, dst, actualAlpha, blocking), this->rows(), this->cols(), Dest::Flags&RowMajorBit);
    }
};

#endif // EIGEN_GENERAL_MATRIX_MATRIX_H
