import {
	Curve,
	Vector3,
	Vector4
} from 'three';
import * as NURBSUtils from '../curves/NURBSUtils.js';

/**
 * NURBS curve object
 *
 * Derives from Curve, overriding getPoint and getTangent.
 *
 * Implementation is based on (x, y [, z=0 [, w=1]]) control points with w=weight.
 *
 **/

class NURBSCurve extends Curve {

	constructor(
		degree,
		knots /* array of reals */,
		controlPoints /* array of Vector(2|3|4) */,
		startKnot /* index in knots */,
		endKnot /* index in knots */
	) {

		super();

		this.degree = degree;
		this.knots = knots;
		this.controlPoints = [];
		// Used by periodic NURBS to remove hidden spans
		this.startKnot = startKnot || 0;
		this.endKnot = endKnot || ( this.knots.length - 1 );

		for ( let i = 0; i < controlPoints.length; ++ i ) {

			// ensure Vector4 for control points
			const point = controlPoints[ i ];
			this.controlPoints[ i ] = new Vector4( point.x, point.y, point.z, point.w );

		}

	}

	getPoint( t, optionalTarget = new Vector3() ) {

		const point = optionalTarget;

		const u = this.knots[ this.startKnot ] + t * ( this.knots[ this.endKnot ] - this.knots[ this.startKnot ] ); // linear mapping t->u

		// following results in (wx, wy, wz, w) homogeneous point
		const hpoint = NURBSUtils.calcBSplinePoint( this.degree, this.knots, this.controlPoints, u );

		if ( hpoint.w !== 1.0 ) {

			// project to 3D space: (wx, wy, wz, w) -> (x, y, z, 1)
			hpoint.divideScalar( hpoint.w );

		}

		return point.set( hpoint.x, hpoint.y, hpoint.z );

	}

	getTangent( t, optionalTarget = new Vector3() ) {

		const tangent = optionalTarget;

		const u = this.knots[ 0 ] + t * ( this.knots[ this.knots.length - 1 ] - this.knots[ 0 ] );
		const ders = NURBSUtils.calcNURBSDerivatives( this.degree, this.knots, this.controlPoints, u, 1 );
		tangent.copy( ders[ 1 ] ).normalize();

		return tangent;

	}

	computeFrenetFrames2( t ){ 
		var C1 = new Vector3(); 
		var C2 = new Vector3(); 

		const tangent = new Vector3(); 
		const normal = new Vector3(); 
		const binormal = new Vector3();
		
		const ders = NURBSUtils.calcNURBSDerivatives( this.degree, this.knots, this.controlPoints, t, 2 );
		C1.copy( ders[ 1 ] ); 
		C2.copy( ders[ 2 ] ); 
		tangent.copy(C1).normalize(); //versore Tangente
		binormal.crossVectors(C1,C2); //versore Binormale 
		binormal.normalize();
		normal.crossVectors(binormal,tangent); //versore normale
		return { tangents: tangent, normals: normal, binormals: binormal }; 
	}

}

export { NURBSCurve };
