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
 *    \file src/code_generation/integrators/irk_lifted_forward_export.cpp
 *    \author Rien Quirynen
 *    \date 2014
 */

#include <acado/code_generation/integrators/irk_export.hpp>
#include <acado/code_generation/integrators/irk_lifted_forward_export.hpp>

using namespace std;

BEGIN_NAMESPACE_ACADO

//
// PUBLIC MEMBER FUNCTIONS:
//

ForwardLiftedIRKExport::ForwardLiftedIRKExport(	UserInteraction* _userInteraction,
									const std::string& _commonHeaderName
									) : ForwardIRKExport( _userInteraction,_commonHeaderName )
{
}

ForwardLiftedIRKExport::ForwardLiftedIRKExport( const ForwardLiftedIRKExport& arg ) : ForwardIRKExport( arg )
{
}


ForwardLiftedIRKExport::~ForwardLiftedIRKExport( )
{
	if ( solver )
		delete solver;
	solver = 0;

	clear( );
}


ForwardLiftedIRKExport& ForwardLiftedIRKExport::operator=( const ForwardLiftedIRKExport& arg ){

    if( this != &arg ){

    	ForwardIRKExport::operator=( arg );
		copy( arg );
    }
    return *this;
}


ExportVariable ForwardLiftedIRKExport::getAuxVariable() const
{
	ExportVariable max;
	if( NX1 > 0 ) {
		max = lin_input.getGlobalExportVariable();
	}
	if( NX2 > 0 || NXA > 0 ) {
		if( rhs.getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = rhs.getGlobalExportVariable();
		}
		if( diffs_rhs.getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = diffs_rhs.getGlobalExportVariable();
		}
		if( diffs_sweep.getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = diffs_sweep.getGlobalExportVariable();
		}
	}
	if( NX3 > 0 ) {
		if( rhs3.getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = rhs3.getGlobalExportVariable();
		}
		if( diffs_rhs3.getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = diffs_rhs3.getGlobalExportVariable();
		}
	}
	uint i;
	for( i = 0; i < outputs.size(); i++ ) {
		if( outputs[i].getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = outputs[i].getGlobalExportVariable();
		}
		if( diffs_outputs[i].getGlobalExportVariable().getDim() >= max.getDim() ) {
			max = diffs_outputs[i].getGlobalExportVariable();
		}
	}
	return max;
}


returnValue ForwardLiftedIRKExport::getDataDeclarations(	ExportStatementBlock& declarations,
												ExportStruct dataStruct
												) const
{
	ForwardIRKExport::getDataDeclarations( declarations, dataStruct );
	
	declarations.addDeclaration( rk_seed,dataStruct );
	declarations.addDeclaration( rk_stageValues,dataStruct );

	declarations.addDeclaration( rk_Xprev,dataStruct );
	declarations.addDeclaration( rk_Uprev,dataStruct );
	declarations.addDeclaration( rk_delta,dataStruct );

    return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::getFunctionDeclarations(	ExportStatementBlock& declarations
													) const
{
	ForwardIRKExport::getFunctionDeclarations( declarations );

    return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::setDifferentialEquation(	const Expression& rhs_ )
{
	int sensGen;
	get( DYNAMIC_SENSITIVITY,sensGen );
	if( rhs_.getDim() > 0 ) {
		OnlineData        dummy0;
		Control           dummy1;
		DifferentialState dummy2;
		AlgebraicState 	  dummy3;
		DifferentialStateDerivative dummy4;
		dummy0.clearStaticCounters();
		dummy1.clearStaticCounters();
		dummy2.clearStaticCounters();
		dummy3.clearStaticCounters();
		dummy4.clearStaticCounters();

		NX2 = rhs_.getDim() - NXA;
		x = DifferentialState("", NX1+NX2, 1);
		z = AlgebraicState("", NXA, 1);
		u = Control("", NU, 1);
		od = OnlineData("", NOD, 1);

		DifferentialEquation f;
		f << rhs_;

		NDX2 = f.getNDX();
		if( NDX2 > 0 && (NDX2 < NX2 || NDX2 > (NX1+NX2)) ) {
			return ACADOERROR( RET_INVALID_OPTION );
		}
		else if( NDX2 > 0 ) NDX2 = NX1+NX2;
		dx = DifferentialStateDerivative("", NDX2, 1);

		DifferentialEquation g;
		for( uint i = 0; i < rhs_.getDim(); i++ ) {
			g << forwardDerivative( rhs_(i), x );
			g << forwardDerivative( rhs_(i), z );
			g << forwardDerivative( rhs_(i), u );
			g << forwardDerivative( rhs_(i), dx );
		}

		DifferentialEquation h;
		if( (ExportSensitivityType)sensGen == INEXACT ) {
			DifferentialState sX("", NX,1), sU("", NU,1);

			// forward sweep
			h << forwardDerivative( rhs_, x, sX ) + forwardDerivative( rhs_, u, sU );
		}

		if( f.getNT() > 0 ) timeDependant = true;

		return (rhs.init( f,"acado_rhs",NX,NXA,NU,NP,NDX,NOD ) &
				diffs_rhs.init( g,"acado_diffs",NX,NXA,NU,NP,NDX,NOD ) &
				diffs_sweep.init( h,"acado_diff_sweep",NX+(NX+NU),NXA,NU,NP,NDX,NOD ) );
	}
	return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::getCode(	ExportStatementBlock& code )
{
	int sensGen;
	get( DYNAMIC_SENSITIVITY, sensGen );
	int mode;
	get( IMPLICIT_INTEGRATOR_MODE, mode );
	int liftMode;
	get( LIFTED_INTEGRATOR_MODE, liftMode );
	if ( (ExportSensitivityType)sensGen != FORWARD && (ExportSensitivityType)sensGen != INEXACT ) ACADOERROR( RET_INVALID_OPTION );
	if( (ImplicitIntegratorMode)mode != LIFTED ) ACADOERROR( RET_INVALID_OPTION );
	if( liftMode > 4 || liftMode == 2 || liftMode < 1 ) ACADOERROR( RET_INVALID_OPTION );
	if( (ExportSensitivityType)sensGen == INEXACT && liftMode != 4 ) ACADOERROR( RET_INVALID_OPTION );

	if( CONTINUOUS_OUTPUT || NX1 > 0 || NX3 > 0 || !equidistantControlGrid() || NDX2 > 0 || NXA > 0 ) ACADOERROR( RET_NOT_IMPLEMENTED_YET );

	int useOMP;
	get(CG_USE_OPENMP, useOMP);
	if ( useOMP ) {
		ExportVariable max = getAuxVariable();
		max.setName( "auxVar" );
		max.setDataStruct( ACADO_LOCAL );
		if( NX2 > 0 || NXA > 0 ) {
			rhs.setGlobalExportVariable( max );
			diffs_rhs.setGlobalExportVariable( max );
			diffs_sweep.setGlobalExportVariable( max );
		}
		if( NX3 > 0 ) {
			rhs3.setGlobalExportVariable( max );
			diffs_rhs3.setGlobalExportVariable( max );
		}
		for( uint i = 0; i < outputs.size(); i++ ) {
			outputs[i].setGlobalExportVariable( max );
			diffs_outputs[i].setGlobalExportVariable( max );
		}

		getDataDeclarations( code, ACADO_LOCAL );

		stringstream s;
		s << "#pragma omp threadprivate( "
				<< max.getFullName() << ", "
				<< rk_ttt.getFullName() << ", "
				<< rk_xxx.getFullName() << ", "
				<< rk_kkk.getFullName() << ", "
				<< rk_diffK.getFullName() << ", "
				<< rk_rhsTemp.getFullName() << ", "
				<< rk_auxSolver.getFullName();
		if( NX1 > 0 ) {
			if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) s << ", " << rk_diffsPrev1.getFullName();
			s << ", " << rk_diffsNew1.getFullName();
		}
		if( NX2 > 0 || NXA > 0 ) {
			s << ", " << rk_A.getFullName();
			if( (ExportSensitivityType)sensGen == INEXACT ) s << ", " << rk_seed.getFullName();
			s << ", " << rk_b.getFullName();
			if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) s << ", " << rk_diffsPrev2.getFullName();
			s << ", " << rk_diffsNew2.getFullName();
			s << ", " << rk_diffsTemp2.getFullName();
			solver->appendVariableNames( s );
		}
		if( NX3 > 0 ) {
			if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) s << ", " << rk_diffsPrev3.getFullName();
			s << ", " << rk_diffsNew3.getFullName();
			s << ", " << rk_diffsTemp3.getFullName();
		}
		s << " )" << endl << endl;
		code.addStatement( s.str().c_str() );
	}

	if( NX1 > 0 ) {
		code.addFunction( lin_input );
		code.addStatement( "\n\n" );
	}
	if( exportRhs ) {
		if( NX2 > 0 || NXA > 0 ) {
			code.addFunction( rhs );
			code.addStatement( "\n\n" );
			code.addFunction( diffs_rhs );
			code.addStatement( "\n\n" );
			code.addFunction( diffs_sweep );
			code.addStatement( "\n\n" );
		}

		if( NX3 > 0 ) {
			code.addFunction( rhs3 );
			code.addStatement( "\n\n" );
			code.addFunction( diffs_rhs3 );
			code.addStatement( "\n\n" );
		}

		if( CONTINUOUS_OUTPUT ) {
			uint i;
			for( i = 0; i < outputs.size(); i++ ) {
				code.addFunction( outputs[i] );
				code.addStatement( "\n\n" );
				code.addFunction( diffs_outputs[i] );
				code.addStatement( "\n\n" );
			}
		}
	}
	if( NX2 > 0 || NXA > 0 ) solver->getCode( code );
	code.addLinebreak(2);

	int measGrid;
	get( MEASUREMENT_GRID, measGrid );

	// export RK scheme
	uint run5;
	std::string tempString;
	
	initializeDDMatrix();
	initializeCoefficients();

	double h = (grid.getLastTime() - grid.getFirstTime())/grid.getNumIntervals();
	DMatrix tmp = AA;
	ExportVariable Ah( "Ah_mat", tmp*=h, STATIC_CONST_REAL );
	code.addDeclaration( Ah );
	code.addLinebreak( 2 );
	// TODO: Ask Milan why this does NOT work properly !!
	Ah = ExportVariable( "Ah_mat", numStages, numStages, STATIC_CONST_REAL, ACADO_LOCAL );

	DVector BB( bb );
	ExportVariable Bh( "Bh_mat", DMatrix( BB*=h ) );

	DVector CC( cc );
	ExportVariable C;
	if( timeDependant ) {
		C = ExportVariable( "C_mat", DMatrix( CC*=(1.0/grid.getNumIntervals()) ), STATIC_CONST_REAL );
		code.addDeclaration( C );
		code.addLinebreak( 2 );
		C = ExportVariable( "C_mat", 1, numStages, STATIC_CONST_REAL, ACADO_LOCAL );
	}

	code.addComment(std::string("Fixed step size:") + toString(h));

	ExportVariable determinant( "det", 1, 1, REAL, ACADO_LOCAL, true );
	integrate.addDeclaration( determinant );

	ExportIndex i( "i" );
	ExportIndex j( "j" );
	ExportIndex k( "k" );
	ExportIndex run( "run" );
	ExportIndex run1( "run1" );
	ExportIndex tmp_index1("tmp_index1");
	ExportIndex tmp_index2("tmp_index2");
	ExportIndex tmp_index3("tmp_index3");
	ExportIndex tmp_index4("tmp_index4");
	ExportIndex k_index("k_index");
	ExportIndex shooting_index("shoot_index");
	ExportVariable tmp_meas("tmp_meas", 1, outputGrids.size(), INT, ACADO_LOCAL);

	ExportVariable numInt( "numInts", 1, 1, INT );
	if( !equidistantControlGrid() ) {
		ExportVariable numStepsV( "numSteps", numSteps, STATIC_CONST_INT );
		code.addDeclaration( numStepsV );
		code.addLinebreak( 2 );
		integrate.addStatement( std::string( "int " ) + numInt.getName() + " = " + numStepsV.getName() + "[" + rk_index.getName() + "];\n" );
	}

	prepareOutputEvaluation( code );

	integrate.addIndex( i );
	integrate.addIndex( j );
	integrate.addIndex( k );
	integrate.addIndex( run );
	integrate.addIndex( run1 );
	integrate.addIndex( tmp_index1 );
	integrate.addIndex( tmp_index2 );
	integrate.addIndex( shooting_index );
	integrate.addIndex( k_index );
	if( rk_outputs.size() > 0 ) integrate.addIndex( tmp_index3 );
	if( rk_outputs.size() > 0 && (grid.getNumIntervals() > 1 || !equidistantControlGrid()) ) {
		integrate.addIndex( tmp_index4 );
	}
	ExportVariable time_tmp( "time_tmp", 1, 1, REAL, ACADO_LOCAL, true );
	if( CONTINUOUS_OUTPUT ) {
		for( run5 = 0; run5 < outputGrids.size(); run5++ ) {
			ExportIndex numMeasTmp( (std::string)"numMeasTmp" + toString(run5) );
			numMeas.push_back( numMeasTmp );
			integrate.addIndex( numMeas[run5] );
		}

		if( (MeasurementGrid)measGrid == ONLINE_GRID ) {
			integrate.addDeclaration( tmp_meas );
			integrate.addDeclaration( polynEvalVar );
			integrate.addDeclaration( time_tmp );
		}

		for( run5 = 0; run5 < outputGrids.size(); run5++ ) {
			integrate.addStatement( numMeas[run5] == 0 );
		}
	}
	integrate << shooting_index.getFullName() << " = " << rk_index.getFullName() << ";\n";
	integrate.addStatement( rk_ttt == DMatrix(grid.getFirstTime()) );
	if( (inputDim-diffsDim) > NX+NXA ) {
		integrate.addStatement( rk_xxx.getCols( NX+NXA,inputDim-diffsDim ) == rk_eta.getCols( NX+NXA+diffsDim,inputDim ) );
		if( (ExportSensitivityType)sensGen == INEXACT ) {
			integrate.addStatement( rk_seed.getCols( NX+NXA+NX+NU,NX+NU+inputDim-diffsDim ) == rk_eta.getCols( NX+NXA+diffsDim,inputDim ) );
		}
	}
	integrate.addLinebreak( );
	if( liftMode == 1 || (liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT) ) {
		integrate.addStatement( rk_delta.getCols( NX,NX+NU ) == rk_eta.getCols( NX+NXA+diffsDim,NX+NXA+diffsDim+NU ) - rk_Uprev.getRow(shooting_index) );
	}

    // integrator loop:
	ExportForLoop tmpLoop( run, 0, grid.getNumIntervals() );
	ExportStatementBlock *loop;
	if( equidistantControlGrid() ) {
		loop = &tmpLoop;
	}
	else {
	    loop = &integrate;
		loop->addStatement( std::string("for(") + run.getName() + " = 0; " + run.getName() + " < " + numInt.getName() + "; " + run.getName() + "++ ) {\n" );
	}

	if( CONTINUOUS_OUTPUT && (MeasurementGrid)measGrid == ONLINE_GRID ) {
		for( run5 = 0; run5 < outputGrids.size(); run5++ ) {
			loop->addStatement( tmp_index1 == numMeas[run5] );
			loop->addStatement( std::string("while( ") + tmp_index1.getName() + " < " + toString(totalMeas[run5]) + " && " + gridVariables[run5].get(0,tmp_index1) + " <= (" + rk_ttt.getFullName() + "+" + toString(1.0/grid.getNumIntervals()) + ") ) {\n" );
			loop->addStatement( tmp_index1 == tmp_index1+1 );
			loop->addStatement( std::string("}\n") );
			loop->addStatement( std::string(tmp_meas.get( 0,run5 )) + " = " + tmp_index1.getName() + " - " + numMeas[run5].getName() + ";\n" );
		}
	}

	if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) {
		// Set rk_diffsPrev:
		loop->addStatement( std::string("if( run > 0 ) {\n") );
		if( NX1 > 0 ) {
			ExportForLoop loopTemp1( i,0,NX1 );
			loopTemp1.addStatement( rk_diffsPrev1.getSubMatrix( i,i+1,0,NX1 ) == rk_eta.getCols( i*NX+NX+NXA,i*NX+NX+NXA+NX1 ) );
			if( NU > 0 ) loopTemp1.addStatement( rk_diffsPrev1.getSubMatrix( i,i+1,NX1,NX1+NU ) == rk_eta.getCols( i*NU+(NX+NXA)*(NX+1),i*NU+(NX+NXA)*(NX+1)+NU ) );
			loop->addStatement( loopTemp1 );
		}
		if( NX2 > 0 ) {
			ExportForLoop loopTemp2( i,0,NX2 );
			loopTemp2.addStatement( rk_diffsPrev2.getSubMatrix( i,i+1,0,NX1+NX2 ) == rk_eta.getCols( i*NX+NX+NXA+NX1*NX,i*NX+NX+NXA+NX1*NX+NX1+NX2 ) );
			if( NU > 0 ) loopTemp2.addStatement( rk_diffsPrev2.getSubMatrix( i,i+1,NX1+NX2,NX1+NX2+NU ) == rk_eta.getCols( i*NU+(NX+NXA)*(NX+1)+NX1*NU,i*NU+(NX+NXA)*(NX+1)+NX1*NU+NU ) );
			loop->addStatement( loopTemp2 );
		}
		if( NX3 > 0 ) {
			ExportForLoop loopTemp3( i,0,NX3 );
			loopTemp3.addStatement( rk_diffsPrev3.getSubMatrix( i,i+1,0,NX ) == rk_eta.getCols( i*NX+NX+NXA+(NX1+NX2)*NX,i*NX+NX+NXA+(NX1+NX2)*NX+NX ) );
			if( NU > 0 ) loopTemp3.addStatement( rk_diffsPrev3.getSubMatrix( i,i+1,NX,NX+NU ) == rk_eta.getCols( i*NU+(NX+NXA)*(NX+1)+(NX1+NX2)*NU,i*NU+(NX+NXA)*(NX+1)+(NX1+NX2)*NU+NU ) );
			loop->addStatement( loopTemp3 );
		}
		loop->addStatement( std::string("}\n") );
	}
	if( liftMode == 1 || (liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT) ) {
		loop->addStatement( rk_delta.getCols( 0,NX ) == rk_eta.getCols( 0,NX ) - rk_Xprev.getRow(shooting_index*grid.getNumIntervals()+run) );
	}

//	// PART 1: The linear input system
//	prepareInputSystem( code );
//	solveInputSystem( loop, i, run1, j, tmp_index1, Ah );

	// PART 2: The fully implicit system
	loop->addStatement( k_index == (shooting_index*grid.getNumIntervals()+run)*(NX+NXA) );
	solveImplicitSystem( loop, i, run1, j, tmp_index1, k_index, Ah, C, determinant, true );

//	// PART 3: The linear output system
//	prepareOutputSystem( code );
//	solveOutputSystem( loop, i, run1, j, tmp_index1, Ah, true );

//	// generate continuous OUTPUT:
//	generateOutput( loop, run, i, tmp_index2, tmp_index3, tmp_meas, time_tmp, NX+NU );

	if( (ExportSensitivityType)sensGen == INEXACT ) {
		evaluateAllStatesImplicitSystem( loop, k_index, Ah, C, run1, j, tmp_index1 );
	}

	// DERIVATIVES wrt the states (IFT):
	if( NX1 > 0 ) {
		ExportForLoop loop4( run1,0,NX1 );
		// PART 1: The linear input system
		sensitivitiesInputSystem( &loop4, run1, i, Bh, true );
		// PART 2: The fully implicit system
		if( (ExportSensitivityType)sensGen != INEXACT ) {
			sensitivitiesImplicitSystem( &loop4, run1, i, j, tmp_index1, tmp_index2, Ah, Bh, determinant, true, 1 );
		}
		else {
			inexactSensitivitiesImplicitSystem( &loop4, run1, i, j, tmp_index1, tmp_index2, k_index, Ah, Bh, determinant, true, 1 );
		}
		// PART 3: The linear output system
		sensitivitiesOutputSystem( &loop4, run1, i, j, k, tmp_index1, tmp_index2, Ah, Bh, true, 1 );
		// generate sensitivities wrt states for continuous output:
		sensitivitiesOutputs( &loop4, run, run1, i, tmp_index1, tmp_index2, tmp_index3, tmp_meas, time_tmp, true, 0 );
		loop->addStatement( loop4 );
	}
	if( NX2 > 0 ) {
		ExportForLoop loop4( run1,NX1,NX1+NX2 );

		// PART 2: The fully implicit system
		if( (ExportSensitivityType)sensGen != INEXACT ) {
			sensitivitiesImplicitSystem( &loop4, run1, i, j, tmp_index1, tmp_index2, Ah, Bh, determinant, true, 2 );
		}
		else {
			inexactSensitivitiesImplicitSystem( &loop4, run1, i, j, tmp_index1, tmp_index2, k_index, Ah, Bh, determinant, true, 2 );
		}
//		// PART 3: The linear output system
//		sensitivitiesOutputSystem( &loop4, run1, i, j, k, tmp_index1, tmp_index2, Ah, Bh, true, 2 );
//		// generate sensitivities wrt states for continuous output:
//		sensitivitiesOutputs( &loop4, run, run1, i, tmp_index1, tmp_index2, tmp_index3, tmp_meas, time_tmp, true, NX1 );

		loop->addStatement( loop4 );
	}
	if( NX3 > 0 ) {
		ExportForLoop loop4( run1,NX1+NX2,NX );
		// PART 3: The linear output system
		sensitivitiesOutputSystem( &loop4, run1, i, j, k, tmp_index1, tmp_index2, Ah, Bh, true, 3 );
		// generate sensitivities wrt states for continuous output:
		sensitivitiesOutputs( &loop4, run, run1, i, tmp_index1, tmp_index2, tmp_index3, tmp_meas, time_tmp, true, NX1+NX2 );
		loop->addStatement( loop4 );
	}


	// DERIVATIVES wrt the control inputs (IFT):
	if( NU > 0 ) {
		ExportForLoop loop5( run1,0,NU );

//		// PART 1: The linear input system
//		sensitivitiesInputSystem( &loop5, run1, i, Bh, false );
		// PART 2: The fully implicit system
		if( (ExportSensitivityType)sensGen != INEXACT ) {
			sensitivitiesImplicitSystem( &loop5, run1, i, j, tmp_index1, tmp_index2, Ah, Bh, determinant, false, 0 );
		}
		else {
			inexactSensitivitiesImplicitSystem( &loop5, run1, i, j, tmp_index1, tmp_index2, k_index, Ah, Bh, determinant, false, 0 );
		}
//		// PART 3: The linear output system
//		sensitivitiesOutputSystem( &loop5, run1, i, j, k, tmp_index1, tmp_index2, Ah, Bh, false, 0 );
//		// generate sensitivities wrt controls for continuous output:
//		sensitivitiesOutputs( &loop5, run, run1, i, tmp_index1, tmp_index2, tmp_index3, tmp_meas, time_tmp, false, 0 );

		loop->addStatement( loop5 );
	}

	if( liftMode == 1 || (liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT) ) {
		loop->addStatement( rk_Xprev.getRow(shooting_index*grid.getNumIntervals()+run) == rk_eta.getCols( 0,NX ) );
	}
	// update rk_eta:
	for( run5 = 0; run5 < NX; run5++ ) {
		loop->addStatement( rk_eta.getCol( run5 ) += rk_kkk.getRow( k_index+run5 )*Bh );
	}
	if( NXA > 0) {
		DMatrix tempCoefs( evaluateDerivedPolynomial( 0.0 ) );
		if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
			loop->addStatement( std::string("if( run == 0 ) {\n") );
		}
		for( run5 = 0; run5 < NXA; run5++ ) {
			loop->addStatement( rk_eta.getCol( NX+run5 ) == rk_kkk.getRow( k_index+NX+run5 )*tempCoefs );
		}
		if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
			loop->addStatement( std::string("}\n") );
		}
	}


	// Computation of the sensitivities using the CHAIN RULE:
	if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) {
		loop->addStatement( std::string( "if( run == 0 ) {\n" ) );
	}
//	// PART 1
//	updateInputSystem(loop, i, j, tmp_index2);

	// PART 2
	updateImplicitSystem(loop, i, j, tmp_index2);

//	// PART 3
//	updateOutputSystem(loop, i, j, tmp_index2);

	if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) {
		loop->addStatement( std::string( "}\n" ) );
		loop->addStatement( std::string( "else {\n" ) );
//		// PART 1
//		propagateInputSystem(loop, i, j, k, tmp_index2);

		// PART 2
		propagateImplicitSystem(loop, i, j, k, tmp_index2);

//		// PART 3
//		propagateOutputSystem(loop, i, j, k, tmp_index2);
	}

	if( rk_outputs.size() > 0 && (grid.getNumIntervals() > 1 || !equidistantControlGrid()) ) {
		propagateOutputs( loop, run, run1, i, j, k, tmp_index1, tmp_index2, tmp_index3, tmp_index4, tmp_meas );
	}

	if( grid.getNumIntervals() > 1 || !equidistantControlGrid() ) {
		loop->addStatement( std::string( "}\n" ) );
	}

//	loop->addStatement( std::string( reset_int.get(0,0) ) + " = 0;\n" );

	for( run5 = 0; run5 < rk_outputs.size(); run5++ ) {
		if( (MeasurementGrid)measGrid == OFFLINE_GRID ) {
			loop->addStatement( numMeas[run5].getName() + " += " + numMeasVariables[run5].get(0,run) + ";\n" );
		}
		else { // ONLINE_GRID
			loop->addStatement( numMeas[run5].getName() + " += " + tmp_meas.get(0,run5) + ";\n" );
		}
	}
	loop->addStatement( rk_ttt += DMatrix(1.0/grid.getNumIntervals()) );

    // end of the integrator loop.
    if( !equidistantControlGrid() ) {
		loop->addStatement( "}\n" );
	}
    else {
    	integrate.addStatement( *loop );
    }
    // PART 1
    if( NX1 > 0 ) {
    	DMatrix zeroR = zeros<double>(1, NX2+NX3);
    	ExportForLoop loop1( i,0,NX1 );
    	loop1.addStatement( rk_eta.getCols( i*NX+NX+NXA+NX1,i*NX+NX+NXA+NX ) == zeroR );
    	integrate.addStatement( loop1 );
    }
    // PART 2
    DMatrix zeroR = zeros<double>(1, NX3);
    if( NX2 > 0 && NX3 > 0 ) {
    	ExportForLoop loop2( i,NX1,NX1+NX2 );
    	loop2.addStatement( rk_eta.getCols( i*NX+NX+NXA+NX1+NX2,i*NX+NX+NXA+NX ) == zeroR );
    	integrate.addStatement( loop2 );
    }
    if( NXA > 0 && NX3 > 0 ) {
    	ExportForLoop loop3( i,NX,NX+NXA );
    	loop3.addStatement( rk_eta.getCols( i*NX+NX+NXA+NX1+NX2,i*NX+NX+NXA+NX ) == zeroR );
    	integrate.addStatement( loop3 );
    }

	if( liftMode == 1 || (liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT) ) {
		integrate.addStatement( rk_Uprev.getRow(shooting_index) == rk_eta.getCols( NX+NXA+diffsDim,NX+NXA+diffsDim+NU ) );
	}

    integrate.addStatement( std::string( "if( " ) + determinant.getFullName() + " < 1e-12 ) {\n" );
    integrate.addStatement( error_code == 2 );
    integrate.addStatement( std::string( "} else if( " ) + determinant.getFullName() + " < 1e-6 ) {\n" );
    integrate.addStatement( error_code == 1 );
    integrate.addStatement( std::string( "} else {\n" ) );
    integrate.addStatement( error_code == 0 );
    integrate.addStatement( std::string( "}\n" ) );

	code.addFunction( integrate );
    code.addLinebreak( 2 );

    return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::solveImplicitSystem( ExportStatementBlock* block, const ExportIndex& index1, const ExportIndex& index2, const ExportIndex& index3, const ExportIndex& tmp_index, const ExportIndex& k_index, const ExportVariable& Ah, const ExportVariable& C, const ExportVariable& det, bool DERIVATIVES )
{
	int sensGen;
	get( DYNAMIC_SENSITIVITY, sensGen);
	if( NX2 > 0 || NXA > 0 ) {
		int liftMode;
		get( LIFTED_INTEGRATOR_MODE, liftMode );

		// FIRST update using term from optimization variables:
		if( liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT ) {
			ExportForLoop loopTemp1( index1,0,NX+NXA );
			loopTemp1.addStatement( index3 == k_index+index1 );
			loopTemp1.addStatement( tmp_index == index3*(NX+NU) );
			ExportForLoop loopTemp2( index2,0,numStages );
			loopTemp2.addStatement( rk_kkk.getElement( index3,index2 ) += rk_delta*rk_diffK.getSubMatrix( tmp_index,tmp_index+NX+NU,index2,index2+1 ) );
			loopTemp1.addStatement( loopTemp2 );
			block->addStatement( loopTemp1 );
		}

		// Perform iteration by system solve:
		if( liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT ) {
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( "if( run == 0 ) {\n" );
			}

			evaluateStatesImplicitSystem( block, k_index, Ah, C, ExportIndex(0), index3, tmp_index );
			block->addFunctionCall( getNameDiffsRHS(), rk_xxx, rk_diffsTemp2 );

			ExportForLoop loop01( index2,0,numStages );
			evaluateInexactMatrix( &loop01, index2, index3, tmp_index, k_index, rk_A, Ah, C, true, DERIVATIVES );
			block->addStatement( loop01 );
//			ExportForLoop loop2( index2,0,numStages*(NX2+NXA) );
//			ExportForLoop loop3( index3,0,numStages*(NX2+NXA) );
//			loop3.addStatement( rk_A_or.getElement(index2,index3) == rk_A.getElement(index2,index3) );
//			loop2.addStatement( loop3 );
//			block->addStatement( loop2 );
			block->addStatement( det.getFullName() + " = " + solver->getNameSolveFunction() + "( " + rk_A.getFullName() + ", " + rk_b.getFullName() + ", " + rk_auxSolver.getFullName() + " );\n" );
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( "}\n else {\n" );
				ExportForLoop loop02( index2,0,numStages );
				evaluateStatesImplicitSystem( &loop02, k_index, Ah, C, index2, index3, tmp_index );
				evaluateRhsImplicitSystem( &loop02, k_index, index2 );
				block->addStatement( loop02 );
				block->addFunctionCall( solver->getNameSolveReuseFunction(),rk_A.getAddress(0,0),rk_b.getAddress(0,0),rk_auxSolver.getAddress(0,0) );
				block->addStatement( "}\n" );
			}
		}
		else if( liftMode == 4 ) {
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( "if( run == 0 ) {\n" );
			}
			ExportForLoop loop01( index2,0,numStages );
			evaluateMatrix( &loop01, index2, index3, tmp_index, k_index, rk_A, Ah, C, true, DERIVATIVES );
			block->addStatement( loop01 );
			block->addStatement( det.getFullName() + " = " + solver->getNameSolveFunction() + "( " + rk_A.getFullName() + ", " + rk_b.getFullName() + ", " + rk_auxSolver.getFullName() + " );\n" );
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( "}\n else {\n" );
				ExportForLoop loop02( index2,0,numStages );
				evaluateStatesImplicitSystem( &loop02, k_index, Ah, C, index2, index3, tmp_index );
				loop02.addFunctionCall( getNameDiffsRHS(), rk_xxx, rk_diffsTemp2.getAddress(index2,0) );
				evaluateRhsImplicitSystem( &loop02, k_index, index2 );
				block->addStatement( loop02 );
				block->addFunctionCall( solver->getNameSolveReuseFunction(),rk_A.getAddress(0,0),rk_b.getAddress(0,0),rk_auxSolver.getAddress(0,0) );
				block->addStatement( "}\n" );
			}
		}
		else {
			ExportForLoop loop1( index2,0,numStages );
			evaluateMatrix( &loop1, index2, index3, tmp_index, k_index, rk_A, Ah, C, true, DERIVATIVES );
			if( liftMode == 1 ) {  // Right-hand side term from update optimization variables:
				for( uint i = 0; i < NX2+NXA; i++ ) {
					loop1.addStatement( rk_b.getRow( index2*(NX2+NXA)+i ) -= rk_diffsTemp2.getSubMatrix( index2,index2+1,i*(NVARS2)+NX1,i*(NVARS2)+NX1+NX2 )*(rk_delta.getCols( NX1,NX1+NX2 ).getTranspose()) );
					loop1.addStatement( rk_b.getRow( index2*(NX2+NXA)+i ) -= rk_diffsTemp2.getSubMatrix( index2,index2+1,i*(NVARS2)+NX1+NX2+NXA,i*(NVARS2)+NX1+NX2+NXA+NU )*(rk_delta.getCols( NX,NX+NU ).getTranspose()) );
				}
			}
			block->addStatement( loop1 );
			block->addStatement( det.getFullName() + " = " + solver->getNameSolveFunction() + "( " + rk_A.getFullName() + ", " + rk_b.getFullName() + ", " + rk_auxSolver.getFullName() + " );\n" );
		}
		ExportForLoop loopTemp( index3,0,numStages );
		loopTemp.addStatement( rk_kkk.getSubMatrix( k_index+NX1,k_index+NX1+NX2,index3,index3+1 ) += rk_b.getRows( index3*NX2,index3*NX2+NX2 ) );											// differential states
		if(NXA > 0) loopTemp.addStatement( rk_kkk.getSubMatrix( k_index+NX,k_index+NX+NXA,index3,index3+1 ) += rk_b.getRows( index3*NXA+numStages*NX2,index3*NXA+numStages*NX2+NXA ) );		// algebraic states
		block->addStatement( loopTemp );

		// IF DEBUG MODE:
		int debugMode;
		get( INTEGRATOR_DEBUG_MODE, debugMode );
		if ( (bool)debugMode == true ) {
			block->addStatement( debug_mat == rk_A );
		}
	}

	return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::evaluateInexactMatrix( ExportStatementBlock* block, const ExportIndex& index1, const ExportIndex& index2, const ExportIndex& tmp_index, const ExportIndex& k_index, const ExportVariable& _rk_A, const ExportVariable& Ah, const ExportVariable& C, bool evaluateB, bool DERIVATIVES )
{
	uint i;

	evaluateStatesImplicitSystem( block, k_index, Ah, C, index1, index2, tmp_index );

	ExportForLoop loop2( index2,0,NX2+NXA );
	loop2.addStatement( tmp_index == index1*(NX2+NXA)+index2 );
	for( i = 0; i < numStages; i++ ) { // differential states
		if( NDX2 == 0 ) {
			loop2.addStatement( _rk_A.getSubMatrix( tmp_index,tmp_index+1,i*NX2,i*NX2+NX2 ) == Ah.getElement( index1,i )*rk_diffsTemp2.getSubMatrix( index2,index2+1,NX1,NX1+NX2 ) );
			loop2.addStatement( std::string( "if( " ) + toString(i) + " == " + index1.getName() + " ) " );
			loop2.addStatement( _rk_A.getElement( tmp_index,index2+i*NX2 ) -= 1 );
		}
		else {
			loop2.addStatement( _rk_A.getSubMatrix( tmp_index,tmp_index+1,i*NX2,i*NX2+NX2 ) == Ah.getElement( index1,i )*rk_diffsTemp2.getSubMatrix( index2,index2+1,NX1,NX1+NX2 ) );
			loop2.addStatement( std::string( "if( " ) + toString(i) + " == " + index1.getName() + " ) {\n" );
			loop2.addStatement( _rk_A.getSubMatrix( tmp_index,tmp_index+1,i*NX2,i*NX2+NX2 ) += rk_diffsTemp2.getSubMatrix( index2,index2+1,NVARS2-NX2,NVARS2 ) );
			loop2.addStatement( std::string( "}\n" ) );
		}
	}
	if( NXA > 0 ) {
		DMatrix zeroM = zeros<double>( 1,NXA );
		for( i = 0; i < numStages; i++ ) { // algebraic states
			loop2.addStatement( std::string( "if( " ) + toString(i) + " == " + index1.getName() + " ) {\n" );
			loop2.addStatement( _rk_A.getSubMatrix( tmp_index,tmp_index+1,numStages*NX2+i*NXA,numStages*NX2+i*NXA+NXA ) == rk_diffsTemp2.getSubMatrix( index2,index2+1,NX1+NX2,NX1+NX2+NXA ) );
			loop2.addStatement( std::string( "}\n else {\n" ) );
			loop2.addStatement( _rk_A.getSubMatrix( tmp_index,tmp_index+1,numStages*NX2+i*NXA,numStages*NX2+i*NXA+NXA ) == zeroM );
			loop2.addStatement( std::string( "}\n" ) );
		}
	}
	block->addStatement( loop2 );
	if( evaluateB ) {
		evaluateRhsImplicitSystem( block, k_index, index1 );
	}

	return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::propagateOutputs(	ExportStatementBlock* block, const ExportIndex& index, const ExportIndex& index0, const ExportIndex& index1,
															const ExportIndex& index2, const ExportIndex& index3, const ExportIndex& tmp_index1, const ExportIndex& tmp_index2,
															const ExportIndex& tmp_index3, const ExportIndex& tmp_index4, const ExportVariable& tmp_meas )
{

	return ACADOERROR( RET_NOT_YET_IMPLEMENTED );
}


returnValue ForwardLiftedIRKExport::prepareInputSystem(	ExportStatementBlock& code )
{

	return ACADOERROR( RET_NOT_YET_IMPLEMENTED );
}


returnValue ForwardLiftedIRKExport::prepareOutputSystem(	ExportStatementBlock& code )
{

	return ACADOERROR( RET_NOT_YET_IMPLEMENTED );
}


returnValue ForwardLiftedIRKExport::sensitivitiesInputSystem( ExportStatementBlock* block, const ExportIndex& index1, const ExportIndex& index2, const ExportVariable& Bh, bool STATES )
{

	return ACADOERROR( RET_NOT_YET_IMPLEMENTED );
}


returnValue ForwardLiftedIRKExport::evaluateAllStatesImplicitSystem( ExportStatementBlock* block, const ExportIndex& k_index, const ExportVariable& Ah, const ExportVariable& C, const ExportIndex& stage, const ExportIndex& i, const ExportIndex& tmp_index )
{
	ExportForLoop loop0( stage,0,numStages );
	ExportForLoop loop1( i, 0, NX1+NX2 );
	loop1.addStatement( rk_stageValues.getCol( stage*(NX+NXA)+i ) == rk_eta.getCol( i ) );
	loop1.addStatement( tmp_index == k_index + i );
	for( uint j = 0; j < numStages; j++ ) {
		loop1.addStatement( rk_stageValues.getCol( stage*(NX+NXA)+i ) += Ah.getElement(stage,j)*rk_kkk.getElement( tmp_index,j ) );
	}
	loop0.addStatement( loop1 );

	ExportForLoop loop3( i, 0, NXA );
	loop3.addStatement( tmp_index == k_index + i + NX );
	loop3.addStatement( rk_stageValues.getCol( stage*(NX+NXA)+NX+i ) == rk_kkk.getElement( tmp_index,stage ) );
	loop0.addStatement( loop3 );
	block->addStatement( loop0 );

	return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::inexactSensitivitiesImplicitSystem( ExportStatementBlock* block, const ExportIndex& index1, const ExportIndex& index2, const ExportIndex& index3, const ExportIndex& tmp_index1, const ExportIndex& tmp_index2, const ExportIndex& k_index, const ExportVariable& Ah, const ExportVariable& Bh, const ExportVariable& det, bool STATES, uint number )
{
	if( NX2 > 0 ) {
		DMatrix zeroM = zeros<double>( NX2+NXA,1 );
		DMatrix tempCoefs( evaluateDerivedPolynomial( 0.0 ) );  // We compute the algebraic variables at the beginning of the shooting interval !
		uint j;

		ExportForLoop loop1( index2,0,numStages );
		loop1.addStatement( rk_seed.getCols(0,NX) == rk_stageValues.getCols(index2*(NX+NXA),index2*(NX+NXA)+NX+NXA) );
		loop1.addStatement( rk_seed.getCols(NX,NX+NX+NU) == zeros<double>(1,NX+NU) );
		if( STATES && number == 1 ) { // TODO: NOT YET IMPLEMENTED (linear input)
//			ExportForLoop loop2( index3,0,NX1 );
//			loop2.addStatement( std::string(rk_rhsTemp.get( index3,0 )) + " = -(" + index3.getName() + " == " + index1.getName() + ");\n" );
//			for( i = 0; i < numStages; i++ ) {
//				loop2.addStatement( rk_rhsTemp.getRow( index3 ) -= rk_diffK.getElement( index3,i )*Ah.getElement(index2,i) );
//			}
//			loop1.addStatement( loop2 );
//			ExportForLoop loop3( index3,0,NX2+NXA );
//			loop3.addStatement( tmp_index1 == index2*(NX2+NXA)+index3 );
//			loop3.addStatement( rk_b.getRow( tmp_index1 ) == rk_diffsTemp2.getSubMatrix( index2,index2+1,index3*(NVARS2),index3*(NVARS2)+NX1 )*rk_rhsTemp.getRows(0,NX1) );
//			if( NDX2 > 0 ) {
//				loop3.addStatement( rk_b.getRow( tmp_index1 ) -= rk_diffsTemp2.getSubMatrix( index2,index2+1,index3*(NVARS2)+NVARS2-NX1-NX2,index3*(NVARS2)+NVARS2-NX2 )*rk_diffK.getSubMatrix( 0,NX1,index2,index2+1 ) );
//			}
//			loop1.addStatement( loop3 );
		}
		else if( STATES && number == 2 ) {
			loop1.addStatement( rk_seed.getCol(NX+index1) -= 1.0 );
			ExportForLoop loop3( index3,0,NX2+NXA );
			loop3.addStatement( tmp_index1 == k_index + index3 );
			loop3.addStatement( tmp_index2 == tmp_index1*(NX+NU) + index1 );
			for( j = 0; j < numStages; j++ ) {
				loop3.addStatement( rk_seed.getCol( NX+index3 ) -= Ah.getElement(index2,j)*rk_diffK.getElement(tmp_index2,j) );
			}
			loop1.addStatement( loop3 );
		}
		else { // ~= STATES
//			ExportForLoop loop2( index3,0,NX1 );
//			loop2.addStatement( rk_rhsTemp.getRow( index3 ) == rk_diffK.getElement( index3,0 )*Ah.getElement(index2,0) );
//			for( i = 1; i < numStages; i++ ) {
//				loop2.addStatement( rk_rhsTemp.getRow( index3 ) += rk_diffK.getElement( index3,i )*Ah.getElement(index2,i) );
//			}
//			loop1.addStatement( loop2 );

			loop1.addStatement( rk_seed.getCol(NX+NX+index1) -= 1.0 );
			ExportForLoop loop3( index3,0,NX2+NXA );
			loop3.addStatement( tmp_index1 == k_index + index3 );
			loop3.addStatement( tmp_index2 == tmp_index1*(NX+NU) + NX + index1 );
			for( j = 0; j < numStages; j++ ) {
				loop3.addStatement( rk_seed.getCol( NX+index3 ) -= Ah.getElement(index2,j)*rk_diffK.getElement(tmp_index2,j) );
			}
			loop1.addStatement( loop3 );
		}
		loop1.addFunctionCall( diffs_sweep.getName(), rk_seed, rk_b.getAddress(index2*(NX2+NXA),0) );
		ExportForLoop loop03( index3,0,NX2+NXA );
		loop03.addStatement( tmp_index1 == k_index + index3 );
		if( STATES ) {
			loop03.addStatement( tmp_index2 == tmp_index1*(NX+NU) + index1 );
		}
		else {
			loop03.addStatement( tmp_index2 == tmp_index1*(NX+NU) + NX + index1 );
		}
		loop03.addStatement( rk_b.getRow( index2*(NX2+NXA)+index3 ) += rk_diffK.getElement(tmp_index2,index2) );
		loop1.addStatement( loop03 );
		block->addStatement( loop1 );

		block->addFunctionCall( solver->getNameSolveReuseFunction(),rk_A.getAddress(0,0),rk_b.getAddress(0,0),rk_auxSolver.getAddress(0,0) );

		// update rk_diffK with the new sensitivities:
		ExportForLoop loop20( index2,0,numStages );
		ExportForLoop loop21( index3,0,NX2 );
		loop21.addStatement( tmp_index1 == (k_index + NX1 + index3)*(NX+NU) );
		if( STATES ) {
			loop21.addStatement( tmp_index2 == tmp_index1+index1 );
		}
		else {
			loop21.addStatement( tmp_index2 == tmp_index1+NX+index1 );
		}
		loop21.addStatement( rk_diffK.getElement(tmp_index2,index2) += rk_b.getRow(index2*NX2+index3) );
		loop20.addStatement( loop21 );
		if( NXA > 0 ) {
			ExportForLoop loop22( index3,0,NXA );
			loop22.addStatement( tmp_index1 == (k_index + NX + index3)*(NX+NU) );
			if( STATES ) {
				loop22.addStatement( tmp_index2 == tmp_index1+index1 );
			}
			else {
				loop22.addStatement( tmp_index2 == tmp_index1+NX+index1 );
			}
			loop22.addStatement( rk_diffK.getElement(tmp_index2,index2) += rk_b.getRow(numStages*NX2+index2*NXA+index3) );
			loop20.addStatement( loop22 );
		}
		block->addStatement( loop20 );

		// update rk_diffsNew with the new sensitivities:
		ExportForLoop loop3( index2,0,NX2 );
		if( STATES && number == 2 ) loop3.addStatement( std::string(rk_diffsNew2.get( index2,index1 )) + " = (" + index2.getName() + " == " + index1.getName() + "-" + toString(NX1) + ");\n" );

		loop3.addStatement( tmp_index1 == (k_index + NX1 + index2)*(NX+NU) );
		loop3.addStatement( tmp_index2 == tmp_index1+index1 );
		if( STATES && number == 2 ) loop3.addStatement( rk_diffsNew2.getElement( index2,index1 ) += rk_diffK.getRow( tmp_index2 )*Bh );
		else if( STATES )	loop3.addStatement( rk_diffsNew2.getElement( index2,index1 ) == rk_diffK.getRow( tmp_index2 )*Bh );
		else		 		loop3.addStatement( rk_diffsNew2.getElement( index2,index1+NX1+NX2 ) == rk_diffK.getRow( tmp_index2+NX )*Bh );
		block->addStatement( loop3 );
		if( NXA > 0 ) {
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( std::string("if( run == 0 ) {\n") );
			}
			ExportForLoop loop4( index2,0,NXA );
			loop4.addStatement( tmp_index1 == (k_index + NX + index2)*(NX+NU) );
			loop4.addStatement( tmp_index2 == tmp_index1+index1 );
			if( STATES ) loop4.addStatement( rk_diffsNew2.getElement( index2+NX2,index1 ) == rk_diffK.getRow( tmp_index2 )*tempCoefs );
			else 		 loop4.addStatement( rk_diffsNew2.getElement( index2+NX2,index1+NX1+NX2 ) == rk_diffK.getRow( tmp_index2+NX )*tempCoefs );
			block->addStatement( loop4 );
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( std::string("}\n") );
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::sensitivitiesImplicitSystem( ExportStatementBlock* block, const ExportIndex& index1, const ExportIndex& index2, const ExportIndex& index3, const ExportIndex& tmp_index1, const ExportIndex& tmp_index2, const ExportVariable& Ah, const ExportVariable& Bh, const ExportVariable& det, bool STATES, uint number )
{
	if( NX2 > 0 ) {
		DMatrix zeroM = zeros<double>( NX2+NXA,1 );
		DMatrix tempCoefs( evaluateDerivedPolynomial( 0.0 ) );  // We compute the algebraic variables at the beginning of the shooting interval !
		uint i;

		ExportForLoop loop1( index2,0,numStages );
		if( STATES && number == 1 ) { // TODO: NOT YET IMPLEMENTED (linear input)
//			ExportForLoop loop2( index3,0,NX1 );
//			loop2.addStatement( std::string(rk_rhsTemp.get( index3,0 )) + " = -(" + index3.getName() + " == " + index1.getName() + ");\n" );
//			for( i = 0; i < numStages; i++ ) {
//				loop2.addStatement( rk_rhsTemp.getRow( index3 ) -= rk_diffK.getElement( index3,i )*Ah.getElement(index2,i) );
//			}
//			loop1.addStatement( loop2 );
//			ExportForLoop loop3( index3,0,NX2+NXA );
//			loop3.addStatement( tmp_index1 == index2*(NX2+NXA)+index3 );
//			loop3.addStatement( rk_b.getRow( tmp_index1 ) == rk_diffsTemp2.getSubMatrix( index2,index2+1,index3*(NVARS2),index3*(NVARS2)+NX1 )*rk_rhsTemp.getRows(0,NX1) );
//			if( NDX2 > 0 ) {
//				loop3.addStatement( rk_b.getRow( tmp_index1 ) -= rk_diffsTemp2.getSubMatrix( index2,index2+1,index3*(NVARS2)+NVARS2-NX1-NX2,index3*(NVARS2)+NVARS2-NX2 )*rk_diffK.getSubMatrix( 0,NX1,index2,index2+1 ) );
//			}
//			loop1.addStatement( loop3 );
		}
		else if( STATES && number == 2 ) {
			for( i = 0; i < NX2+NXA; i++ ) {
				loop1.addStatement( rk_b.getRow( index2*(NX2+NXA)+i ) == zeroM.getRow( 0 ) - rk_diffsTemp2.getElement( index2,index1+i*(NVARS2) ) );
			}
		}
		else { // ~= STATES
//			ExportForLoop loop2( index3,0,NX1 );
//			loop2.addStatement( rk_rhsTemp.getRow( index3 ) == rk_diffK.getElement( index3,0 )*Ah.getElement(index2,0) );
//			for( i = 1; i < numStages; i++ ) {
//				loop2.addStatement( rk_rhsTemp.getRow( index3 ) += rk_diffK.getElement( index3,i )*Ah.getElement(index2,i) );
//			}
//			loop1.addStatement( loop2 );
			ExportForLoop loop3( index3,0,NX2+NXA );
			loop3.addStatement( tmp_index1 == index2*(NX2+NXA)+index3 );
			loop3.addStatement( tmp_index2 == index1+index3*(NVARS2) );
			loop3.addStatement( rk_b.getRow( tmp_index1 ) == zeroM.getRow( 0 ) - rk_diffsTemp2.getElement( index2,tmp_index2+NX1+NX2+NXA ) );
//			loop3.addStatement( rk_b.getRow( tmp_index1 ) -= rk_diffsTemp2.getSubMatrix( index2,index2+1,index3*(NVARS2),index3*(NVARS2)+NX1 )*rk_rhsTemp.getRows(0,NX1) );
//			if( NDX2 > 0 ) {
//				loop3.addStatement( rk_b.getRow( tmp_index1 ) -= rk_diffsTemp2.getSubMatrix( index2,index2+1,index3*(NVARS2)+NVARS2-NX1-NX2,index3*(NVARS2)+NVARS2-NX2 )*rk_diffK.getSubMatrix( 0,NX1,index2,index2+1 ) );
//			}
			loop1.addStatement( loop3 );
		}
		block->addStatement( loop1 );

		block->addFunctionCall( solver->getNameSolveReuseFunction(),rk_A.getAddress(0,0),rk_b.getAddress(0,0),rk_auxSolver.getAddress(0,0) );

		// update rk_diffK with the new sensitivities:
		ExportForLoop loop2( index2,0,numStages );
		loop2.addStatement( rk_diffK.getSubMatrix(NX1,NX1+NX2,index2,index2+1) == rk_b.getRows(index2*NX2,index2*NX2+NX2) );
		loop2.addStatement( rk_diffK.getSubMatrix(NX,NX+NXA,index2,index2+1) == rk_b.getRows(numStages*NX2+index2*NXA,numStages*NX2+index2*NXA+NXA) );
		block->addStatement( loop2 );
		// update rk_diffsNew with the new sensitivities:
		ExportForLoop loop3( index2,0,NX2 );
		if( STATES && number == 2 ) loop3.addStatement( std::string(rk_diffsNew2.get( index2,index1 )) + " = (" + index2.getName() + " == " + index1.getName() + "-" + toString(NX1) + ");\n" );

		if( STATES && number == 2 ) loop3.addStatement( rk_diffsNew2.getElement( index2,index1 ) += rk_diffK.getRow( NX1+index2 )*Bh );
		else if( STATES )	loop3.addStatement( rk_diffsNew2.getElement( index2,index1 ) == rk_diffK.getRow( NX1+index2 )*Bh );
		else		 		loop3.addStatement( rk_diffsNew2.getElement( index2,index1+NX1+NX2 ) == rk_diffK.getRow( NX1+index2 )*Bh );
		block->addStatement( loop3 );
		if( NXA > 0 ) {
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( std::string("if( run == 0 ) {\n") );
			}
			ExportForLoop loop4( index2,0,NXA );
			if( STATES ) loop4.addStatement( rk_diffsNew2.getElement( index2+NX2,index1 ) == rk_diffK.getRow( NX+index2 )*tempCoefs );
			else 		 loop4.addStatement( rk_diffsNew2.getElement( index2+NX2,index1+NX1+NX2 ) == rk_diffK.getRow( NX+index2 )*tempCoefs );
			block->addStatement( loop4 );
			if( !equidistantControlGrid() || grid.getNumIntervals() > 1 ) {
				block->addStatement( std::string("}\n") );
			}
		}
	}

	return SUCCESSFUL_RETURN;
}


returnValue ForwardLiftedIRKExport::sensitivitiesOutputSystem( ExportStatementBlock* block, const ExportIndex& index1, const ExportIndex& index2, const ExportIndex& index3, const ExportIndex& index4, const ExportIndex& tmp_index1, const ExportIndex& tmp_index2, const ExportVariable& Ah, const ExportVariable& Bh, bool STATES, uint number )
{

	return ACADOERROR( RET_NOT_YET_IMPLEMENTED );
}


returnValue ForwardLiftedIRKExport::sensitivitiesOutputs( ExportStatementBlock* block, const ExportIndex& index0,
		const ExportIndex& index1, const ExportIndex& index2, const ExportIndex& tmp_index1, const ExportIndex& tmp_index2,
		const ExportIndex& tmp_index3, const ExportVariable& tmp_meas, const ExportVariable& time_tmp, bool STATES, uint base )
{

	return ACADOERROR( RET_NOT_YET_IMPLEMENTED );
}


returnValue ForwardLiftedIRKExport::setup( )
{
	ForwardIRKExport::setup();

	integrate = ExportFunction( "integrate", rk_eta );
	uint i;
	for( i = 0; i < rk_outputs.size(); i++ ) {
		integrate.addArgument( rk_outputs[i] );
	}
	integrate.addArgument( rk_index );
	integrate.setReturnValue( error_code );
	integrate.doc( "Performs the integration and sensitivity propagation for one shooting interval." );
	integrate.addLinebreak( );	// TO MAKE SURE IT GETS EXPORTED

	int useOMP;
	get(CG_USE_OPENMP, useOMP);
	ExportStruct structWspace;
	structWspace = useOMP ? ACADO_LOCAL : ACADO_WORKSPACE;

	uint timeDep = 0;
	if( timeDependant ) timeDep = 1;

	int sensGen;
	get( DYNAMIC_SENSITIVITY, sensGen);
	if( (ExportSensitivityType)sensGen == INEXACT ) {
		rk_seed = ExportVariable( "rk_seed", 1, NX+NX+NU+NXA+NU+NOD+NDX+timeDep, REAL, structWspace );
		rk_diffsTemp2 = ExportVariable( "rk_diffsTemp2", NX2+NXA, NVARS2, REAL, structWspace );
		rk_stageValues = ExportVariable( "rk_stageValues", 1, numStages*(NX+NXA), REAL, structWspace );
	}

	rk_kkk = ExportVariable( "rk_Ktraj", N*grid.getNumIntervals()*(NX+NXA), numStages, REAL, ACADO_VARIABLES );
	if( (ExportSensitivityType)sensGen == INEXACT ) {
		rk_diffK = ExportVariable( "rk_diffKtraj", N*grid.getNumIntervals()*(NX+NXA)*(NX+NU), numStages, REAL, ACADO_VARIABLES );
	}

	int liftMode;
	get( LIFTED_INTEGRATOR_MODE, liftMode );
	if( liftMode == 1 || (liftMode == 4 && (ExportSensitivityType)sensGen == INEXACT) ) {
		rk_Xprev = ExportVariable( "rk_Xprev", N*grid.getNumIntervals(), NX, REAL, ACADO_VARIABLES );
		rk_Uprev = ExportVariable( "rk_Uprev", N, NU, REAL, ACADO_VARIABLES );
		rk_delta = ExportVariable( "rk_delta", 1, NX+NU, REAL, ACADO_WORKSPACE );
	}

    return SUCCESSFUL_RETURN;
}



// PROTECTED:


CLOSE_NAMESPACE_ACADO

// end of file.