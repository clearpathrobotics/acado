/*
 *    This file is part of ACADO Toolkit.
 *
 *    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
 *    Copyright (C) 2008-2014 by Boris Houska, Hans Joachim Ferreau,
 *    Milan Vukov, Rien Quirynen, KU Leuven.
 *    Developed within the Optimization in Engineering Center (OPTEC)
 *    under supervision of Moritz Diehl. All rights reserved.
 *
 *    ACADO Toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADO Toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *    \file src/nlp_solver/scp_merit_function.cpp
 *    \author Boris Houska, Hans Joachim Ferreau
 *
 */


#include <acado/nlp_solver/scp_merit_function.hpp>



BEGIN_NAMESPACE_ACADO


//
// PUBLIC MEMBER FUNCTIONS:
//

SCPmeritFunction::SCPmeritFunction( ) : AlgorithmicBase( )
{
}


SCPmeritFunction::SCPmeritFunction( UserInteraction* _userInteraction ) : AlgorithmicBase( _userInteraction )
{
}


SCPmeritFunction::SCPmeritFunction( const SCPmeritFunction& rhs ) : AlgorithmicBase( rhs )
{
}


SCPmeritFunction::~SCPmeritFunction( )
{
}


SCPmeritFunction& SCPmeritFunction::operator=( const SCPmeritFunction& rhs )
{
    if ( this != &rhs )
    {
		AlgorithmicBase::operator=( rhs );
    }

    return *this;
}


SCPmeritFunction* SCPmeritFunction::clone( ) const
{
	return new SCPmeritFunction( *this );
}


returnValue SCPmeritFunction::evaluate(	double alpha,
										const OCPiterate& iter,
										BandedCP& cp,
										SCPevaluation& eval,
										double& result
										)
{
    if( fabs( alpha ) >= EPS )
	{
		result = INFTY;

        eval.clearDynamicDiscretization( );

		OCPiterate iterTest = iter;
		iterTest.applyStep( cp.deltaX,alpha );

		if ( eval.evaluate( iterTest,cp ) != SUCCESSFUL_RETURN )
			return SUCCESSFUL_RETURN;
    }

    result = eval.getObjectiveValue( );


    const double kappa = 1.2;
    DMatrix tmp;

//     acadoPrintf("result1 = %.16e \n", result );

//     acadoPrintf("cp.dynResiduum = \n");
//     cp.dynResiduum.print();

//     acadoPrintf("cp.lambdaDynamic = \n");
//     cp.lambdaDynamic.print();

    if( eval.isDynamicNLP( ) == BT_TRUE ){
        //(cp.lambdaDynamic.getAbsolute()^cp.dynResiduum.getAbsolute()).getSubBlock( 0, 0, tmp, 1, 1 );
		cp.lambdaDynamic.getAbsolute(cp_lambdaDynamic_getAbsolute);
		cp.dynResiduum.getAbsolute(cp_dynResiduum_getAbsolute);
		BlockMatrix::matT_mat_mul(cp_lambdaDynamic_getAbsolute, cp_dynResiduum_getAbsolute, temporary0);
		temporary0.getSubBlock( 0, 0, tmp, 1, 1 );
        result += kappa*tmp(0,0);
    }

//     acadoPrintf("result2 = %.16e \n", result );
//
//    --------
//
//
//      acadoPrintf("cp.lambdaBound = \n");
//     (cp.lambdaBound.absolute()       ).print();
//      acadoPrintf("upperBoundRes = \n");
//     (cp.upperBoundResiduum.negative()).print();
//      acadoPrintf("lowerBoundRes = \n");
//     (cp.lowerBoundResiduum.negative()).print();

	
	cp.lambdaBound.getAbsolute(cp_lambdaBound_getAbsolute);
	cp.upperBoundResiduum.getNegative(cp_upperBoundResiduum_getNegative);
	BlockMatrix::matT_mat_mul(cp_lambdaBound_getAbsolute, cp_upperBoundResiduum_getNegative, temporary1);
	temporary1.getSubBlock( 0, 0, tmp, 1, 1 );
    result -= kappa*tmp(0,0);

//     acadoPrintf("result3 = %.16e \n", result );

	cp.lowerBoundResiduum.getPositive(cp_lowerBoundResiduum_getPositive);
	BlockMatrix::matT_mat_mul(cp_lambdaBound_getAbsolute, cp_lowerBoundResiduum_getPositive, temporary2);
	temporary2.getSubBlock( 0, 0, tmp, 1, 1 );
	result += kappa*tmp(0,0);

//     acadoPrintf("result4 = %.16e \n", result );

    // --------

	cp.lambdaConstraint.getAbsolute(cp_lambdaConstraint_getAbsolute);
	cp.upperConstraintResiduum.getNegative(cp_upperConstraintResiduum_getNegative);
	BlockMatrix::matT_mat_mul(cp_lambdaConstraint_getAbsolute, cp_upperConstraintResiduum_getNegative, temporary3);
	temporary3.getSubBlock( 0, 0, tmp, 1, 1 );
    result -= kappa*tmp(0,0);

//      acadoPrintf("result5 = %.16e \n", result );
//
//      acadoPrintf("cp.lambdaConstraint = \n");
//      cp.lambdaConstraint.print();
//
//      acadoPrintf("cp.lambdaConstraint = \n");
//      cp.lowerConstraintResiduum.print();

	cp.lowerConstraintResiduum.getPositive(cp_lowerConstraintResiduum_getPositive);
	BlockMatrix::matT_mat_mul(cp_lambdaConstraint_getAbsolute, cp_lowerConstraintResiduum_getPositive, temporary4);
	temporary4.getSubBlock( 0, 0, tmp, 1, 1 );

	result += kappa*tmp(0,0);

//     acadoPrintf("result6 = %.16e \n", result );

    return SUCCESSFUL_RETURN;
}




CLOSE_NAMESPACE_ACADO

// end of file.
