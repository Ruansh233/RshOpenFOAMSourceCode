/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  dev                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General m4 macros

changecom(//)changequote([,]) dnl>
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// geometric parameters defined by the user

scale   1;

// input parameter
define(X, 0.05)                
define(Y, 0.05)        
define(Z, 0.2)                

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// parametric description of the geometry

vertices
(
    //buttom
    (-X  -Y  -Z) vlabel(q01)
    ( X  -Y  -Z) vlabel(q02)
    ( X   Y  -Z) vlabel(q03)
    (-X   Y  -Z) vlabel(q04)

    (-X  -Y   Z) vlabel(q05)
    ( X  -Y   Z) vlabel(q06)
    ( X   Y   Z) vlabel(q07)
    (-X   Y   Z) vlabel(q08) 

);


blocks
( 
    // x y z of the block must equal to the coordinate of points
    hex ( q01 q02 q03 q04 q05 q06 q07 q08) (5 5 20) simpleGrading (1 1 1) 
);


edges 
(

);  

boundary
(

   wall_minX
   {
      type patch;
      faces
      (
      (q01 q04 q08 q05) 
      );
   }

   wall_maxX
   {
      type patch;
      faces
      (
      (q02 q03 q07 q06) 
      );
   }

   wall_minY
   {
      type patch;
      faces
      (
      (q01 q02 q06 q05) 
      );
   }

   wall_maxY
   {
      type patch;
      faces
      (
      (q04 q03 q07 q08) 
      );
   }

   inlet_minZ
   {
      type patch;
      faces    
      (
      (q01 q02 q03 q04) 
      );
   }

   outlet_maxZ
   {
      type patch;
      faces    
      (
      (q05 q06 q07 q08) 
      );
   }

);

// ************************************************************************* //

