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

scale 0.001;

// input parameter
define(a, 10)
define(L, 20)                                

// point coordinate
define(x1, calc(-1/2*a))
define(y1, calc(-1/2*a))

define(x2, calc(1/2*a))
define(y2, calc(-1/2*a))

define(x3, calc(1/2*a))
define(y3, calc(1/2*a))

define(x4, calc(-1/2*a))
define(y4, calc(1/2*a))

// Number of nodes
define(Nsq, 10);
define(NL,  20);                                          

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// parametric description of the geometry

vertices
(
    //in
    (x1  y1  0) vlabel(q01)
    (x2  y2  0) vlabel(q02)

    (x3  y3  0) vlabel(q03)
    (x4  y4  0) vlabel(q04)

    //out
    (x1  y1  L) vlabel(q05)
    (x2  y2  L) vlabel(q06) 

    (x3  y3  L) vlabel(q07)
    (x4  y4  L) vlabel(q08)

);


blocks
( 
    // x y z of the block must equal to the coordinate of points
    hex (q01 q02 q03 q04 q05 q06 q07 q08) (Nsq Nsq NL) simpleGrading (1 1 1)
);


// edges 
// (
// );  

boundary
(
    block1_in
    {
       type patch;
       faces    
       (
        (q01 q02 q03 q04) 
       );
    }

    block1_out
    {
       type patch;
       faces    
       (
        (q05 q06 q07 q08) 
       );
    }

    block1_wall1
    {
       type wall;
       faces    
       (
        (q01 q05 q08 q04) 
       );
    }

    block1_wall2
    {
       type wall;
       faces    
       (
        (q04 q08 q07 q03) 
       );
    }

    block1_wall3
    {
       type wall;
       faces    
       (
        (q03 q07 q06 q02) 
       );
    }

    block1_wall4
    {
       type wall;
       faces    
       (
        (q01 q05 q06 q02) 
       );
    }


);

// ************************************************************************* //

