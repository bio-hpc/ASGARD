// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2007-2010 Benoit Jacob <jacob.benoit.1@gmail.com>
// Copyright (C) 2008-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
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

#ifndef EIGEN_DENSEBASE_H
#define EIGEN_DENSEBASE_H

/** \class DenseBase
  * \ingroup Core_Module
  *
  * \brief Base class for all dense matrices, vectors, and arrays
  *
  * This class is the base that is inherited by all dense objects (matrix, vector, arrays,
  * and related expression types). The common Eigen API for dense objects is contained in this class.
  *
  * \param Derived is the derived type, e.g., a matrix type or an expression.
  *
  * \sa \ref TopicClassHierarchy
  */
template<typename Derived> class DenseBase
#ifndef EIGEN_PARSED_BY_DOXYGEN
  : public ei_special_scalar_op_base<Derived,typename ei_traits<Derived>::Scalar,
                                     typename NumTraits<typename ei_traits<Derived>::Scalar>::Real>
#else
  : public DenseCoeffsBase<Derived>
#endif // not EIGEN_PARSED_BY_DOXYGEN
{
  public:
    using ei_special_scalar_op_base<Derived,typename ei_traits<Derived>::Scalar,
                typename NumTraits<typename ei_traits<Derived>::Scalar>::Real>::operator*;

    class InnerIterator;

    typedef typename ei_traits<Derived>::StorageKind StorageKind;
    typedef typename ei_traits<Derived>::Index Index; /**< The type of indices */
    typedef typename ei_traits<Derived>::Scalar Scalar;
    typedef typename ei_packet_traits<Scalar>::type PacketScalar;
    typedef typename NumTraits<Scalar>::Real RealScalar;

    typedef DenseCoeffsBase<Derived> Base;
    using Base::derived;
    using Base::const_cast_derived;
    using Base::rows;
    using Base::cols;
    using Base::size;
    using Base::rowIndexByOuterInner;
    using Base::colIndexByOuterInner;
    using Base::coeff;
    using Base::coeffByOuterInner;
    using Base::packet;
    using Base::packetByOuterInner;
    using Base::writePacket;
    using Base::writePacketByOuterInner;
    using Base::coeffRef;
    using Base::coeffRefByOuterInner;
    using Base::copyCoeff;
    using Base::copyCoeffByOuterInner;
    using Base::copyPacket;
    using Base::copyPacketByOuterInner;
    using Base::operator();
    using Base::operator[];
    using Base::x;
    using Base::y;
    using Base::z;
    using Base::w;
    using Base::stride;
    using Base::innerStride;
    using Base::outerStride;
    using Base::rowStride;
    using Base::colStride;
    typedef typename Base::CoeffReturnType CoeffReturnType;

    enum {

      RowsAtCompileTime = ei_traits<Derived>::RowsAtCompileTime,
        /**< The number of rows at compile-time. This is just a copy of the value provided
          * by the \a Derived type. If a value is not known at compile-time,
          * it is set to the \a Dynamic constant.
          * \sa MatrixBase::rows(), MatrixBase::cols(), ColsAtCompileTime, SizeAtCompileTime */

      ColsAtCompileTime = ei_traits<Derived>::ColsAtCompileTime,
        /**< The number of columns at compile-time. This is just a copy of the value provided
          * by the \a Derived type. If a value is not known at compile-time,
          * it is set to the \a Dynamic constant.
          * \sa MatrixBase::rows(), MatrixBase::cols(), RowsAtCompileTime, SizeAtCompileTime */


      SizeAtCompileTime = (ei_size_at_compile_time<ei_traits<Derived>::RowsAtCompileTime,
                                                   ei_traits<Derived>::ColsAtCompileTime>::ret),
        /**< This is equal to the number of coefficients, i.e. the number of
          * rows times the number of columns, or to \a Dynamic if this is not
          * known at compile-time. \sa RowsAtCompileTime, ColsAtCompileTime */

      MaxRowsAtCompileTime = ei_traits<Derived>::MaxRowsAtCompileTime,
        /**< This value is equal to the maximum possible number of rows that this expression
          * might have. If this expression might have an arbitrarily high number of rows,
          * this value is set to \a Dynamic.
          *
          * This value is useful to know when evaluating an expression, in order to determine
          * whether it is possible to avoid doing a dynamic memory allocation.
          *
          * \sa RowsAtCompileTime, MaxColsAtCompileTime, MaxSizeAtCompileTime
          */

      MaxColsAtCompileTime = ei_traits<Derived>::MaxColsAtCompileTime,
        /**< This value is equal to the maximum possible number of columns that this expression
          * might have. If this expression might have an arbitrarily high number of columns,
          * this value is set to \a Dynamic.
          *
          * This value is useful to know when evaluating an expression, in order to determine
          * whether it is possible to avoid doing a dynamic memory allocation.
          *
          * \sa ColsAtCompileTime, MaxRowsAtCompileTime, MaxSizeAtCompileTime
          */

      MaxSizeAtCompileTime = (ei_size_at_compile_time<ei_traits<Derived>::MaxRowsAtCompileTime,
                                                      ei_traits<Derived>::MaxColsAtCompileTime>::ret),
        /**< This value is equal to the maximum possible number of coefficients that this expression
          * might have. If this expression might have an arbitrarily high number of coefficients,
          * this value is set to \a Dynamic.
          *
          * This value is useful to know when evaluating an expression, in order to determine
          * whether it is possible to avoid doing a dynamic memory allocation.
          *
          * \sa SizeAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime
          */

      IsVectorAtCompileTime = ei_traits<Derived>::MaxRowsAtCompileTime == 1
                           || ei_traits<Derived>::MaxColsAtCompileTime == 1,
        /**< This is set to true if either the number of rows or the number of
          * columns is known at compile-time to be equal to 1. Indeed, in that case,
          * we are dealing with a column-vector (if there is only one column) or with
          * a row-vector (if there is only one row). */

      Flags = ei_traits<Derived>::Flags,
        /**< This stores expression \ref flags flags which may or may not be inherited by new expressions
          * constructed from this one. See the \ref flags "list of flags".
          */

      IsRowMajor = int(Flags) & RowMajorBit, /**< True if this expression has row-major storage order. */

      InnerSizeAtCompileTime = int(IsVectorAtCompileTime) ? SizeAtCompileTime
                             : int(IsRowMajor) ? ColsAtCompileTime : RowsAtCompileTime,

      CoeffReadCost = ei_traits<Derived>::CoeffReadCost,
        /**< This is a rough measure of how expensive it is to read one coefficient from
          * this expression.
          */

      InnerStrideAtCompileTime = ei_inner_stride_at_compile_time<Derived>::ret,
      OuterStrideAtCompileTime = ei_outer_stride_at_compile_time<Derived>::ret
    };

    /** \returns the number of nonzero coefficients which is in practice the number
      * of stored coefficients. */
    inline Index nonZeros() const { return size(); }
    /** \returns true if either the number of rows or the number of columns is equal to 1.
      * In other words, this function returns
      * \code rows()==1 || cols()==1 \endcode
      * \sa rows(), cols(), IsVectorAtCompileTime. */

    /** \returns the outer size.
      *
      * \note For a vector, this returns just 1. For a matrix (non-vector), this is the major dimension
      * with respect to the storage order, i.e., the number of columns for a column-major matrix,
      * and the number of rows for a row-major matrix. */
    Index outerSize() const
    {
      return IsVectorAtCompileTime ? 1
           : int(IsRowMajor) ? this->rows() : this->cols();
    }

    /** \returns the inner size.
      *
      * \note For a vector, this is just the size. For a matrix (non-vector), this is the minor dimension
      * with respect to the storage order, i.e., the number of rows for a column-major matrix,
      * and the number of columns for a row-major matrix. */
    Index innerSize() const
    {
      return IsVectorAtCompileTime ? this->size()
           : int(IsRowMajor) ? this->cols() : this->rows();
    }

    /** Only plain matrices/arrays, not expressions, may be resized; therefore the only useful resize methods are
      * Matrix::resize() and Array::resize(). The present method only asserts that the new size equals the old size, and does
      * nothing else.
      */
    void resize(Index size)
    {
      EIGEN_ONLY_USED_FOR_DEBUG(size);
      ei_assert(size == this->size()
                && "DenseBase::resize() does not actually allow to resize.");
    }
    /** Only plain matrices/arrays, not expressions, may be resized; therefore the only useful resize methods are
      * Matrix::resize() and Array::resize(). The present method only asserts that the new size equals the old size, and does
      * nothing else.
      */
    void resize(Index rows, Index cols)
    {
      EIGEN_ONLY_USED_FOR_DEBUG(rows);
      EIGEN_ONLY_USED_FOR_DEBUG(cols);
      ei_assert(rows == this->rows() && cols == this->cols()
                && "DenseBase::resize() does not actually allow to resize.");
    }

#ifndef EIGEN_PARSED_BY_DOXYGEN

    /** \internal Represents a matrix with all coefficients equal to one another*/
    typedef CwiseNullaryOp<ei_scalar_constant_op<Scalar>,Derived> ConstantReturnType;
    /** \internal Represents a vector with linearly spaced coefficients that allows sequential access only. */
    typedef CwiseNullaryOp<ei_linspaced_op<Scalar,false>,Derived> SequentialLinSpacedReturnType;
    /** \internal Represents a vector with linearly spaced coefficients that allows random access. */
    typedef CwiseNullaryOp<ei_linspaced_op<Scalar,true>,Derived> RandomAccessLinSpacedReturnType;
    /** \internal the return type of MatrixBase::eigenvalues() */
    typedef Matrix<typename NumTraits<typename ei_traits<Derived>::Scalar>::Real, ei_traits<Derived>::ColsAtCompileTime, 1> EigenvaluesReturnType;

#endif // not EIGEN_PARSED_BY_DOXYGEN

    /** Copies \a other into *this. \returns a reference to *this. */
    template<typename OtherDerived>
    Derived& operator=(const DenseBase<OtherDerived>& other);

    /** Special case of the template operator=, in order to prevent the compiler
      * from generating a default operator= (issue hit with g++ 4.1)
      */
    Derived& operator=(const DenseBase& other);

    template<typename OtherDerived>
    Derived& operator=(const EigenBase<OtherDerived> &other);

    template<typename OtherDerived>
    Derived& operator+=(const EigenBase<OtherDerived> &other);

    template<typename OtherDerived>
    Derived& operator-=(const EigenBase<OtherDerived> &other);

    template<typename OtherDerived>
    Derived& operator=(const ReturnByValue<OtherDerived>& func);

#ifndef EIGEN_PARSED_BY_DOXYGEN
    /** Copies \a other into *this without evaluating other. \returns a reference to *this. */
    template<typename OtherDerived>
    Derived& lazyAssign(const DenseBase<OtherDerived>& other);
#endif // not EIGEN_PARSED_BY_DOXYGEN

    CommaInitializer<Derived> operator<< (const Scalar& s);

    template<unsigned int Added,unsigned int Removed>
    const Flagged<Derived, Added, Removed> flagged() const;

    template<typename OtherDerived>
    CommaInitializer<Derived> operator<< (const DenseBase<OtherDerived>& other);

    Eigen::Transpose<Derived> transpose();
    const Eigen::Transpose<Derived> transpose() const;
    void transposeInPlace();
#ifndef EIGEN_NO_DEBUG
  protected:
    template<typename OtherDerived>
    void checkTransposeAliasing(const OtherDerived& other) const;
  public:
#endif

    VectorBlock<Derived> segment(Index start, Index size);
    const VectorBlock<Derived> segment(Index start, Index size) const;

    VectorBlock<Derived> head(Index size);
    const VectorBlock<Derived> head(Index size) const;

    VectorBlock<Derived> tail(Index size);
    const VectorBlock<Derived> tail(Index size) const;

    template<int Size> VectorBlock<Derived,Size> head(void);
    template<int Size> const VectorBlock<Derived,Size> head() const;

    template<int Size> VectorBlock<Derived,Size> tail();
    template<int Size> const VectorBlock<Derived,Size> tail() const;

    template<int Size> VectorBlock<Derived,Size> segment(Index start);
    template<int Size> const VectorBlock<Derived,Size> segment(Index start) const;

    Diagonal<Derived,0> diagonal();
    const Diagonal<Derived,0> diagonal() const;

    template<int Index> Diagonal<Derived,Index> diagonal();
    template<int Index> const Diagonal<Derived,Index> diagonal() const;

    Diagonal<Derived, Dynamic> diagonal(Index index);
    const Diagonal<Derived, Dynamic> diagonal(Index index) const;

    template<unsigned int Mode> TriangularView<Derived, Mode> part();
    template<unsigned int Mode> const TriangularView<Derived, Mode> part() const;

    template<unsigned int Mode> TriangularView<Derived, Mode> triangularView();
    template<unsigned int Mode> const TriangularView<Derived, Mode> triangularView() const;

    template<unsigned int UpLo> SelfAdjointView<Derived, UpLo> selfadjointView();
    template<unsigned int UpLo> const SelfAdjointView<Derived, UpLo> selfadjointView() const;

    static const ConstantReturnType
    Constant(Index rows, Index cols, const Scalar& value);
    static const ConstantReturnType
    Constant(Index size, const Scalar& value);
    static const ConstantReturnType
    Constant(const Scalar& value);

    static const SequentialLinSpacedReturnType
    LinSpaced(Sequential_t, Index size, const Scalar& low, const Scalar& high);
    static const RandomAccessLinSpacedReturnType
    LinSpaced(Index size, const Scalar& low, const Scalar& high);
    static const SequentialLinSpacedReturnType
    LinSpaced(Sequential_t, const Scalar& low, const Scalar& high);
    static const RandomAccessLinSpacedReturnType
    LinSpaced(const Scalar& low, const Scalar& high);

    template<typename CustomNullaryOp>
    static const CwiseNullaryOp<CustomNullaryOp, Derived>
    NullaryExpr(Index rows, Index cols, const CustomNullaryOp& func);
    template<typename CustomNullaryOp>
    static const CwiseNullaryOp<CustomNullaryOp, Derived>
    NullaryExpr(Index size, const CustomNullaryOp& func);
    template<typename CustomNullaryOp>
    static const CwiseNullaryOp<CustomNullaryOp, Derived>
    NullaryExpr(const CustomNullaryOp& func);

    static const ConstantReturnType Zero(Index rows, Index cols);
    static const ConstantReturnType Zero(Index size);
    static const ConstantReturnType Zero();
    static const ConstantReturnType Ones(Index rows, Index cols);
    static const ConstantReturnType Ones(Index size);
    static const ConstantReturnType Ones();

    void fill(const Scalar& value);
    Derived& setConstant(const Scalar& value);
    Derived& setLinSpaced(Index size, const Scalar& low, const Scalar& high);
    Derived& setLinSpaced(const Scalar& low, const Scalar& high);
    Derived& setZero();
    Derived& setOnes();
    Derived& setRandom();

    template<typename OtherDerived>
    bool isApprox(const DenseBase<OtherDerived>& other,
                  RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;
    bool isMuchSmallerThan(const RealScalar& other,
                           RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;
    template<typename OtherDerived>
    bool isMuchSmallerThan(const DenseBase<OtherDerived>& other,
                           RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;

    bool isApproxToConstant(const Scalar& value, RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;
    bool isConstant(const Scalar& value, RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;
    bool isZero(RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;
    bool isOnes(RealScalar prec = NumTraits<Scalar>::dummy_precision()) const;

    inline Derived& operator*=(const Scalar& other);
    inline Derived& operator/=(const Scalar& other);

    /** \returns the matrix or vector obtained by evaluating this expression.
      *
      * Notice that in the case of a plain matrix or vector (not an expression) this function just returns
      * a const reference, in order to avoid a useless copy.
      */
    EIGEN_STRONG_INLINE typename ei_eval<Derived>::type eval() const
    {
      // Even though MSVC does not honor strong inlining when the return type
      // is a dynamic matrix, we desperately need strong inlining for fixed
      // size types on MSVC.
      return typename ei_eval<Derived>::type(derived());
    }

    template<typename OtherDerived>
    void swap(DenseBase<OtherDerived> EIGEN_REF_TO_TEMPORARY other);

    inline const NestByValue<Derived> nestByValue() const;
    inline const ForceAlignedAccess<Derived> forceAlignedAccess() const;
    inline ForceAlignedAccess<Derived> forceAlignedAccess();
    template<bool Enable> inline const typename ei_meta_if<Enable,ForceAlignedAccess<Derived>,Derived&>::ret forceAlignedAccessIf() const;
    template<bool Enable> inline typename ei_meta_if<Enable,ForceAlignedAccess<Derived>,Derived&>::ret forceAlignedAccessIf();

    Scalar sum() const;
    Scalar mean() const;
    Scalar trace() const;

    Scalar prod() const;

    typename ei_traits<Derived>::Scalar minCoeff() const;
    typename ei_traits<Derived>::Scalar maxCoeff() const;

    typename ei_traits<Derived>::Scalar minCoeff(Index* row, Index* col) const;
    typename ei_traits<Derived>::Scalar maxCoeff(Index* row, Index* col) const;

    typename ei_traits<Derived>::Scalar minCoeff(Index* index) const;
    typename ei_traits<Derived>::Scalar maxCoeff(Index* index) const;

    template<typename BinaryOp>
    typename ei_result_of<BinaryOp(typename ei_traits<Derived>::Scalar)>::type
    redux(const BinaryOp& func) const;

    template<typename Visitor>
    void visit(Visitor& func) const;

    inline const WithFormat<Derived> format(const IOFormat& fmt) const;

/////////// Array module ///////////

    bool all(void) const;
    bool any(void) const;
    Index count() const;

    const VectorwiseOp<Derived,Horizontal> rowwise() const;
    VectorwiseOp<Derived,Horizontal> rowwise();
    const VectorwiseOp<Derived,Vertical> colwise() const;
    VectorwiseOp<Derived,Vertical> colwise();

    static const CwiseNullaryOp<ei_scalar_random_op<Scalar>,Derived> Random(Index rows, Index cols);
    static const CwiseNullaryOp<ei_scalar_random_op<Scalar>,Derived> Random(Index size);
    static const CwiseNullaryOp<ei_scalar_random_op<Scalar>,Derived> Random();

    template<typename ThenDerived,typename ElseDerived>
    const Select<Derived,ThenDerived,ElseDerived>
    select(const DenseBase<ThenDerived>& thenMatrix,
           const DenseBase<ElseDerived>& elseMatrix) const;

    template<typename ThenDerived>
    inline const Select<Derived,ThenDerived, typename ThenDerived::ConstantReturnType>
    select(const DenseBase<ThenDerived>& thenMatrix, typename ThenDerived::Scalar elseScalar) const;

    template<typename ElseDerived>
    inline const Select<Derived, typename ElseDerived::ConstantReturnType, ElseDerived >
    select(typename ElseDerived::Scalar thenScalar, const DenseBase<ElseDerived>& elseMatrix) const;

    template<int p> RealScalar lpNorm() const;

    template<int RowFactor, int ColFactor>
    const Replicate<Derived,RowFactor,ColFactor> replicate() const;
    const Replicate<Derived,Dynamic,Dynamic> replicate(Index rowFacor,Index colFactor) const;

    Eigen::Reverse<Derived, BothDirections> reverse();
    const Eigen::Reverse<Derived, BothDirections> reverse() const;
    void reverseInPlace();

#define EIGEN_CURRENT_STORAGE_BASE_CLASS Eigen::DenseBase
#   include "../plugins/BlockMethods.h"
#   ifdef EIGEN_DENSEBASE_PLUGIN
#     include EIGEN_DENSEBASE_PLUGIN
#   endif
#undef EIGEN_CURRENT_STORAGE_BASE_CLASS

#ifdef EIGEN2_SUPPORT

    Block<Derived> corner(CornerType type, Index cRows, Index cCols);
    const Block<Derived> corner(CornerType type, Index cRows, Index cCols) const;
    template<int CRows, int CCols>
    Block<Derived, CRows, CCols> corner(CornerType type);
    template<int CRows, int CCols>
    const Block<Derived, CRows, CCols> corner(CornerType type) const;

#endif // EIGEN2_SUPPORT


    // disable the use of evalTo for dense objects with a nice compilation error
    template<typename Dest> inline void evalTo(Dest& ) const
    {
      EIGEN_STATIC_ASSERT((ei_is_same_type<Dest,void>::ret),THE_EVAL_EVALTO_FUNCTION_SHOULD_NEVER_BE_CALLED_FOR_DENSE_OBJECTS);
    }

  protected:
    /** Default constructor. Do nothing. */
    DenseBase()
    {
      /* Just checks for self-consistency of the flags.
       * Only do it when debugging Eigen, as this borders on paranoiac and could slow compilation down
       */
#ifdef EIGEN_INTERNAL_DEBUGGING
      EIGEN_STATIC_ASSERT(ei_are_flags_consistent<Flags>::ret,
                          INVALID_MATRIXBASE_TEMPLATE_PARAMETERS)
      EIGEN_STATIC_ASSERT((EIGEN_IMPLIES(MaxRowsAtCompileTime==1 && MaxColsAtCompileTime!=1, int(IsRowMajor))
                        && EIGEN_IMPLIES(MaxColsAtCompileTime==1 && MaxRowsAtCompileTime!=1, int(!IsRowMajor))),
                          INVALID_STORAGE_ORDER_FOR_THIS_VECTOR_EXPRESSION)
#endif
    }

  private:
    explicit DenseBase(int);
    DenseBase(int,int);
    template<typename OtherDerived> explicit DenseBase(const DenseBase<OtherDerived>&);
};

#endif // EIGEN_DENSEBASE_H
