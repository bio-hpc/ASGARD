// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2006-2008 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2009-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#ifndef EIGEN_TRANSPOSE_H
#define EIGEN_TRANSPOSE_H

/** \class Transpose
  * \ingroup Core_Module
  *
  * \brief Expression of the transpose of a matrix
  *
  * \param MatrixType the type of the object of which we are taking the transpose
  *
  * This class represents an expression of the transpose of a matrix.
  * It is the return type of MatrixBase::transpose() and MatrixBase::adjoint()
  * and most of the time this is the only way it is used.
  *
  * \sa MatrixBase::transpose(), MatrixBase::adjoint()
  */
template<typename MatrixType>
struct ei_traits<Transpose<MatrixType> > : ei_traits<MatrixType>
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename ei_nested<MatrixType>::type MatrixTypeNested;
  typedef typename ei_unref<MatrixTypeNested>::type _MatrixTypeNested;
  typedef typename ei_traits<MatrixType>::StorageKind StorageKind;
  typedef typename ei_traits<MatrixType>::XprKind XprKind;
  enum {
    RowsAtCompileTime = MatrixType::ColsAtCompileTime,
    ColsAtCompileTime = MatrixType::RowsAtCompileTime,
    MaxRowsAtCompileTime = MatrixType::MaxColsAtCompileTime,
    MaxColsAtCompileTime = MatrixType::MaxRowsAtCompileTime,
    Flags = int(_MatrixTypeNested::Flags & ~NestByRefBit) ^ RowMajorBit,
    CoeffReadCost = _MatrixTypeNested::CoeffReadCost,
    InnerStrideAtCompileTime = ei_inner_stride_at_compile_time<MatrixType>::ret,
    OuterStrideAtCompileTime = ei_outer_stride_at_compile_time<MatrixType>::ret
  };
};

template<typename MatrixType, typename StorageKind> class TransposeImpl;

template<typename MatrixType> class Transpose
  : public TransposeImpl<MatrixType,typename ei_traits<MatrixType>::StorageKind>
{
  public:

    typedef typename TransposeImpl<MatrixType,typename ei_traits<MatrixType>::StorageKind>::Base Base;
    EIGEN_GENERIC_PUBLIC_INTERFACE(Transpose)

    inline Transpose(const MatrixType& matrix) : m_matrix(matrix) {}

    EIGEN_INHERIT_ASSIGNMENT_OPERATORS(Transpose)

    inline Index rows() const { return m_matrix.cols(); }
    inline Index cols() const { return m_matrix.rows(); }

    /** \returns the nested expression */
    const typename ei_cleantype<typename MatrixType::Nested>::type&
    nestedExpression() const { return m_matrix; }

    /** \returns the nested expression */
    typename ei_cleantype<typename MatrixType::Nested>::type&
    nestedExpression() { return m_matrix.const_cast_derived(); }

  protected:
    typename MatrixType::Nested m_matrix;
};

template<typename MatrixType, bool HasDirectAccess = ei_has_direct_access<MatrixType>::ret>
struct ei_TransposeImpl_base
{
  typedef typename ei_dense_xpr_base<Transpose<MatrixType> >::type type;
};

template<typename MatrixType>
struct ei_TransposeImpl_base<MatrixType, false>
{
  typedef typename ei_dense_xpr_base<Transpose<MatrixType> >::type type;
};

template<typename MatrixType> class TransposeImpl<MatrixType,Dense>
  : public ei_TransposeImpl_base<MatrixType>::type
{
  public:

    typedef typename ei_TransposeImpl_base<MatrixType>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(Transpose<MatrixType>)

    inline Index innerStride() const { return derived().nestedExpression().innerStride(); }
    inline Index outerStride() const { return derived().nestedExpression().outerStride(); }
    inline Scalar* data() { return derived().nestedExpression().data(); }
    inline const Scalar* data() const { return derived().nestedExpression().data(); }

    inline Scalar& coeffRef(Index row, Index col)
    {
      return const_cast_derived().nestedExpression().coeffRef(col, row);
    }

    inline Scalar& coeffRef(Index index)
    {
      return const_cast_derived().nestedExpression().coeffRef(index);
    }

    inline CoeffReturnType coeff(Index row, Index col) const
    {
      return derived().nestedExpression().coeff(col, row);
    }

    inline CoeffReturnType coeff(Index index) const
    {
      return derived().nestedExpression().coeff(index);
    }

    template<int LoadMode>
    inline const PacketScalar packet(Index row, Index col) const
    {
      return derived().nestedExpression().template packet<LoadMode>(col, row);
    }

    template<int LoadMode>
    inline void writePacket(Index row, Index col, const PacketScalar& x)
    {
      const_cast_derived().nestedExpression().template writePacket<LoadMode>(col, row, x);
    }

    template<int LoadMode>
    inline const PacketScalar packet(Index index) const
    {
      return derived().nestedExpression().template packet<LoadMode>(index);
    }

    template<int LoadMode>
    inline void writePacket(Index index, const PacketScalar& x)
    {
      const_cast_derived().nestedExpression().template writePacket<LoadMode>(index, x);
    }
};

/** \returns an expression of the transpose of *this.
  *
  * Example: \include MatrixBase_transpose.cpp
  * Output: \verbinclude MatrixBase_transpose.out
  *
  * \warning If you want to replace a matrix by its own transpose, do \b NOT do this:
  * \code
  * m = m.transpose(); // bug!!! caused by aliasing effect
  * \endcode
  * Instead, use the transposeInPlace() method:
  * \code
  * m.transposeInPlace();
  * \endcode
  * which gives Eigen good opportunities for optimization, or alternatively you can also do:
  * \code
  * m = m.transpose().eval();
  * \endcode
  *
  * \sa transposeInPlace(), adjoint() */
template<typename Derived>
inline Transpose<Derived>
DenseBase<Derived>::transpose()
{
  return derived();
}

/** This is the const version of transpose().
  *
  * Make sure you read the warning for transpose() !
  *
  * \sa transposeInPlace(), adjoint() */
template<typename Derived>
inline const Transpose<Derived>
DenseBase<Derived>::transpose() const
{
  return derived();
}

/** \returns an expression of the adjoint (i.e. conjugate transpose) of *this.
  *
  * Example: \include MatrixBase_adjoint.cpp
  * Output: \verbinclude MatrixBase_adjoint.out
  *
  * \warning If you want to replace a matrix by its own adjoint, do \b NOT do this:
  * \code
  * m = m.adjoint(); // bug!!! caused by aliasing effect
  * \endcode
  * Instead, use the adjointInPlace() method:
  * \code
  * m.adjointInPlace();
  * \endcode
  * which gives Eigen good opportunities for optimization, or alternatively you can also do:
  * \code
  * m = m.adjoint().eval();
  * \endcode
  *
  * \sa adjointInPlace(), transpose(), conjugate(), class Transpose, class ei_scalar_conjugate_op */
template<typename Derived>
inline const typename MatrixBase<Derived>::AdjointReturnType
MatrixBase<Derived>::adjoint() const
{
  return this->transpose();
}

/***************************************************************************
* "in place" transpose implementation
***************************************************************************/

template<typename MatrixType,
  bool IsSquare = (MatrixType::RowsAtCompileTime == MatrixType::ColsAtCompileTime) && MatrixType::RowsAtCompileTime!=Dynamic>
struct ei_inplace_transpose_selector;

template<typename MatrixType>
struct ei_inplace_transpose_selector<MatrixType,true> { // square matrix
  static void run(MatrixType& m) {
    m.template triangularView<StrictlyUpper>().swap(m.transpose());
  }
};

template<typename MatrixType>
struct ei_inplace_transpose_selector<MatrixType,false> { // non square matrix
  static void run(MatrixType& m) {
    if (m.rows()==m.cols())
      m.template triangularView<StrictlyUpper>().swap(m.transpose());
    else
      m = m.transpose().eval();
  }
};

/** This is the "in place" version of transpose(): it replaces \c *this by its own transpose.
  * Thus, doing
  * \code
  * m.transposeInPlace();
  * \endcode
  * has the same effect on m as doing
  * \code
  * m = m.transpose().eval();
  * \endcode
  * and is faster and also safer because in the latter line of code, forgetting the eval() results
  * in a bug caused by aliasing.
  *
  * Notice however that this method is only useful if you want to replace a matrix by its own transpose.
  * If you just need the transpose of a matrix, use transpose().
  *
  * \note if the matrix is not square, then \c *this must be a resizable matrix.
  *
  * \sa transpose(), adjoint(), adjointInPlace() */
template<typename Derived>
inline void DenseBase<Derived>::transposeInPlace()
{
  ei_inplace_transpose_selector<Derived>::run(derived());
}

/***************************************************************************
* "in place" adjoint implementation
***************************************************************************/

/** This is the "in place" version of adjoint(): it replaces \c *this by its own transpose.
  * Thus, doing
  * \code
  * m.adjointInPlace();
  * \endcode
  * has the same effect on m as doing
  * \code
  * m = m.adjoint().eval();
  * \endcode
  * and is faster and also safer because in the latter line of code, forgetting the eval() results
  * in a bug caused by aliasing.
  *
  * Notice however that this method is only useful if you want to replace a matrix by its own adjoint.
  * If you just need the adjoint of a matrix, use adjoint().
  *
  * \note if the matrix is not square, then \c *this must be a resizable matrix.
  *
  * \sa transpose(), adjoint(), transposeInPlace() */
template<typename Derived>
inline void MatrixBase<Derived>::adjointInPlace()
{
  derived() = adjoint().eval();
}

#ifndef EIGEN_NO_DEBUG

// The following is to detect aliasing problems in most common cases.

template<typename BinOp,typename NestedXpr,typename Rhs>
struct ei_blas_traits<SelfCwiseBinaryOp<BinOp,NestedXpr,Rhs> >
 : ei_blas_traits<NestedXpr>
{
  typedef SelfCwiseBinaryOp<BinOp,NestedXpr,Rhs> XprType;
  static inline const XprType extract(const XprType& x) { return x; }
};

template<bool DestIsTransposed, typename OtherDerived>
struct ei_check_transpose_aliasing_compile_time_selector
{
  enum { ret = ei_blas_traits<OtherDerived>::IsTransposed != DestIsTransposed
  };
};

template<bool DestIsTransposed, typename BinOp, typename DerivedA, typename DerivedB>
struct ei_check_transpose_aliasing_compile_time_selector<DestIsTransposed,CwiseBinaryOp<BinOp,DerivedA,DerivedB> >
{
  enum { ret =    ei_blas_traits<DerivedA>::IsTransposed != DestIsTransposed
               || ei_blas_traits<DerivedB>::IsTransposed != DestIsTransposed
  };
};

template<typename Scalar, bool DestIsTransposed, typename OtherDerived>
struct ei_check_transpose_aliasing_run_time_selector
{
  static bool run(const Scalar* dest, const OtherDerived& src)
  {
    return (ei_blas_traits<OtherDerived>::IsTransposed != DestIsTransposed) && (dest!=0 && dest==(Scalar*)ei_extract_data(src));
  }
};

template<typename Scalar, bool DestIsTransposed, typename BinOp, typename DerivedA, typename DerivedB>
struct ei_check_transpose_aliasing_run_time_selector<Scalar,DestIsTransposed,CwiseBinaryOp<BinOp,DerivedA,DerivedB> >
{
  static bool run(const Scalar* dest, const CwiseBinaryOp<BinOp,DerivedA,DerivedB>& src)
  {
    return ((ei_blas_traits<DerivedA>::IsTransposed != DestIsTransposed) && (dest!=0 && dest==(Scalar*)ei_extract_data(src.lhs())))
        || ((ei_blas_traits<DerivedB>::IsTransposed != DestIsTransposed) && (dest!=0 && dest==(Scalar*)ei_extract_data(src.rhs())));
  }
};

// the following selector, checkTransposeAliasing_impl, based on MightHaveTransposeAliasing,
// is because when the condition controlling the assert is known at compile time, ICC emits a warning.
// This is actually a good warning: in expressions that don't have any transposing, the condition is
// known at compile time to be false, and using that, we can avoid generating the code of the assert again
// and again for all these expressions that don't need it.

template<typename Derived, typename OtherDerived,
         bool MightHaveTransposeAliasing
                 = ei_check_transpose_aliasing_compile_time_selector
                     <ei_blas_traits<Derived>::IsTransposed,OtherDerived>::ret
        >
struct checkTransposeAliasing_impl
{
    static void run(const Derived& dst, const OtherDerived& other)
    {
        ei_assert((!ei_check_transpose_aliasing_run_time_selector
                      <typename Derived::Scalar,ei_blas_traits<Derived>::IsTransposed,OtherDerived>
                      ::run(ei_extract_data(dst), other))
          && "aliasing detected during tranposition, use transposeInPlace() "
             "or evaluate the rhs into a temporary using .eval()");

    }
};

template<typename Derived, typename OtherDerived>
struct checkTransposeAliasing_impl<Derived, OtherDerived, false>
{
    static void run(const Derived&, const OtherDerived&)
    {
    }
};


template<typename Derived>
template<typename OtherDerived>
void DenseBase<Derived>::checkTransposeAliasing(const OtherDerived& other) const
{
    checkTransposeAliasing_impl<Derived, OtherDerived>::run(derived(), other);
}
#endif

#endif // EIGEN_TRANSPOSE_H
