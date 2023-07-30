
/*
//
// peptide.scad
// Copyright (c) 2019 Robert T. Miller.  All rights reserved.
// This file is part of the Biopython distribution and governed by your
// choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
// Please see the LICENSE file that should have been included as part of this
// package.
//
// This is the support file to build an OpenSCAD (http://www.openscad.org/) model
// of a protein from internal coordinates.  The resulting model may be constructed
// on a 3D printer.
//
// data matrices should be appended below to form a program ready to load into
// the OpenSCAD application.
//
//  The protein_scale value used throughout is the second element of the
//    protein[] array appended below.
//    This is the value supplied when generating the data for build units per
//    PDB angstrom.
//    You may wish to modify it here to adjust the appearance of the model in
//    terms of atom sphere or bond cylinder diameter, however the bond lengths
//    are fixed with the supplied value when the data matrices are generated.
//    Atom sphere and bond cylinder radii may be individually adjusted below as
//    well.
//
//  $fn (fragment number) is an OpenSCAD parameter controlling the smoothness
//    of the model surface.  Smaller values will render faster, but yield more
//    'blocky' models.
//
//  This is intended to be a working example, you are encouraged to modify the
//    OpenSCAD subroutines below to generate a model to your liking.  For more
//    information, start with http://www.openscad.org/cheatsheet/index.html
//
//  Note especially the hedronDispatch() subroutine below: here you may select
//    hedra based on residue, sequence position, and class (hedron atoms) for
//    special handling.  Also see the per hedron render options in the hedra[]
//    array.
//
//  If you modify this file, you may find it useful to generate the data
//    matrices without this OpenSCAD code by calling write_SCAD() with the
//    includeCode=False option, then use the OpenSCAD 'include <>' facility at
//    the end of your modified OpenSCAD program.
*/

rotate([-90,0,0])  // convenient for default location (no N-Ca-C start coordinates)
    chain(protein);   // this is the main subroutine call to  build the structure

// top-level OpenSCAD $fn for visible surfaces.  Rotatable bonds use $fn=8
// inside, regardless of this setting.
$fn = 0;  // 0 yields OpenSCAD default of 30.  $n=8 should print with minimal support

tubes=false;     // style: render atoms and bonds as constant diameter cylinders, preferred for rotatable bonds / h-bonds
support=false;   // enable print-in-place internal support for rotatable bonds
// N.B. rotatable bonds must be parallel to build plate for internal support
// structures to be generated correctly by slicer

// output parameters
atomScale=1.0;  // 0.8 better for rotatable bonds
defaultAtomRadius = 0.77;  // used if tubes = true

bondRadius = (tubes ? defaultAtomRadius * atomScale : 0.4);
jBondRadius = defaultAtomRadius * atomScale;  // radius for rotatable bonds

// general printer, slicer, print settings
layerHeight=0.15;  // must match slicer setting for print-in-place support
clearance=0.3;     // sliding clearance - can be smaller (0.2) if not doing print-in-place
pClearance=0.2;    // press-fit clearance (magnets for h-bonds)
shim=0.05;         // extra to make OpenSCAD surfaces distinct in difference()
nozzleDiameter=0.4;

// need one magnet for each side of hydrogen bond, suggest 3mm x 5mm e.g. from eBay
// use compass to identify poles if you care, North pointing (red) repels compass North pointing
magR=3/2;    // magnet radius
magL=5;      // magnet length

// for $fn=8 which works nice on fdm printer
oRot = 22.5;              // 45/2, rotate to make fn=8 spheres and cylinders flat on build plate
apmFac = cos(180/8);      // apothem factor - multiply by radius for center to octagon side distance
octSide = 2* tan(180/8);  // multiply by radius to get length of octagon side
// for values of $fn:
fnRot = ($fn ? 90-(180/$fn) : 90-(180/30));

bondLenFac = 0.6;         // fraction of bond length to extend from atom for each arm of hedron in join

hblen = 1.97;             // hydrogen bond length

wall = 3*nozzleDiameter;
joinerStep = 1;           // radius difference between rotatable bond axle and end knob inside bond cylinder

caTop = false;     // only make top of N_C-alpha_C hedron plus C-beta (see hedron() and hedron_dispatch() examples)

/*
//
// Generate a sphere to represent an atom.
// Colour and size determined for the atom covalent radius specified by the
//   parameter 'a' by lookup in the atomData table below, then scaled by the
//   supplied parameter 'scal'.
//
// scal : protein_scale
// clr : additional radius if used to create clearance for rotatable bonds
//
*/
module atom(a,scal,clr=0)
{
    ad = atomData[search([a],atomData)[0]];
    color(ad[1]) {
        rotate([0,0,fnRot]) sphere(r=((ad[2]*atomScale)*scal)+clr);
    }
}

/*
//
// a hedron (below) may be 'reversed' in terms of the order of its two bonds;
// this function fixes the ordering
//
*/
function hFlip(h,rev) =
        //   yes reversed                                     :  not reversed
        //    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
        //  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
    (rev ? [ h[2], h[0], h[5], h[3], h[8], h[6], h[10], h[9] ] : [ h[0], h[2], h[3], h[5], h[6], h[8],  h[9], h[10] ]);
    // h[1] = angle2 for both cases


/*
//
// generate the male or female interior cylinders of a rotating bond
//
*/
module joinUnit(cOuterLen, cOuterRad, cInnerLen, cInnerRad, male=false) {
    if (male) {
        rotate([0,0,oRot]) {
            cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=8);
            cylinder(h=cOuterLen, r=cOuterRad, center=false, $fn=8);
        }
    } else {
        rotate([0,0,fnRot]) {
            cylinder(h=cInnerLen, r=cInnerRad, center=false, $fn=30);
            cylinder(h=cOuterLen, r=cOuterRad, center=false, $fn=30);
        }
    }
}

/*
//
// create a rotatable bond
//
// supportSel : 0 for no support, 1 or 2 for support on top or bottom (needed
// for reversed hedra)
//
*/
module joiner(bondlen, scal, male=0, ver=0, supportSelect=0) {  // ver = differentiate joiner part lengths to guide assembly, but not used
    lenfac = bondLenFac;
    jClr = clearance+0.05;

    cOuterRad = (jBondRadius * scal) - (2*wall + (male ? jClr/2 : -jClr/2));
    cInnerRad = cOuterRad - joinerStep;  // m/f jClr already in cOuterRad;  - (male ? 0 : -0*jClr/2);

    hArmLen = (bondlen * lenfac);
    lenClr = 0.6*jClr;  // length clearance applied to male and female both, so effective clearance is 2x this value
    cOuterLen = hArmLen * lenfac + (ver ? 0.5 : - 0.5) - (wall+ (male ? lenClr*2 : -lenClr*2  ));

    joinerOffset = (hArmLen * (1 - lenfac)) + (male ? lenClr : -lenClr) - (ver ? 1 : 0);

    i=supportSelect-1;
    oside = cOuterRad*octSide;
    wid = oside+2*wall+4*jClr+1;

    if (male) {
        rotate([0,180,0])
        translate([0,0,-(bondlen-joinerOffset)]) {
            difference() {
                joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad, male=true);
                if (supportSelect) {
                    rotate([0,0,i*180]) {
                        translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                                cube([oside+2*shim,layerHeight+shim,cOuterLen+2*shim],center=true);
                        }
                    }
                }
            }
            if (supportSelect) {
                rotate([0,0,i*180]) {
                    translate([0,(cOuterRad*apmFac)-0.5*layerHeight,cOuterLen/2]) {
                        for (j=[0:1]) {
                            rotate([0,(j?60:-60),0])
                                cube([wid,layerHeight,2*nozzleDiameter],center=true);
                        }
                    }
                }
            }
        }
    } else {
        translate([0,0,joinerOffset]) {
            joinUnit(cOuterLen, cOuterRad, bondlen, cInnerRad);
            if (supportSelect) {  // extra gap top and bottom because filament sags
                supHeight = max(5*layerHeight,2*(cOuterRad-cOuterRad*apmFac));  // double because center=true below
                for(j=[0:1]) {
                    rotate([0,0,j*180]) {
                        translate([0,(cOuterRad*apmFac),cOuterLen/2]) {
                            cube([oside+2*shim,supHeight+shim,cOuterLen+2*shim],center=true);
                        }
                    }
                }
            }
        }
    }
}


/*
//
// create bond with different options (regular, skinny, h-bond atom, rotatable
// male or female
//
//  parameters:
//  bl : bond length
//  br : bond radius
//  scal : protein_scale
//  key : option symbols defined below
//  atm : atomic element symbol, used for color and radius by atom() routine above
//  ver : make rotatable bonds slightly different based on value; currently unused
//  supporSel : enable print-in-place support for rotatable bonds
//
*/

// option symbols - these names generated in BioPython code so avoid changing without thought
StdBond = 1;
FemaleJoinBond = 2;
MaleJoinBond = 3;
SkinnyBond = 4;        // Calpha - Cbeta bond cylinder needs to be skinny for clearance with rotating bonds
HBond = 5;             // make room inside atom/bond to insert magnet to appropriate depth

module bond(bl, br, scal, key, atm, ver, supportSel=0) {

    br = (key == FemaleJoinBond ? jBondRadius * scal : br)  * (key == SkinnyBond ? 0.65 : 1);   // bond radius smaller for skinnyBond
    bl = (key == FemaleJoinBond ? bl * bondLenFac : bl);  // make female joiner shorter
    if (key == MaleJoinBond) { // male join is direct solid, others need difference()
        joiner(bl, scal, male = true, ver = ver, supportSelect=supportSel);
    } else {  // regular bond / skinny / h-bond / female join
        bhblen = bl +(hblen/2 * scal);
        rotate([0,0,fnRot]) {
            difference() {
                union() {
                    cylinder(h=bl,r=br,center=false);
                    if (key == HBond) {  // make extension collar for h-bond magnet
                        rotate([0,0,oRot-fnRot]) cylinder(h=bhblen-1,r=(magR + clearance +wall),center=false, $fn=8);
                    }
                }
                atom(atm,scal,-clearance);  // remove overlap with atom to clear area for female join
                if (key == HBond) {     // make space to insert magnet inside bond cylinder
                    translate([0,0,(bhblen-magL)-pClearance])
                        cylinder(h=magL+pClearance+shim, r=magR+pClearance, center=false, $fn=8);
                }
            }
        }
    }
}

/*
//
// Generate a 'hedron', one plane of 3 points, consisting of 3 atoms joined by
//   two bonds.
//   Defined as bond length - bond angle - bond length
//
// In some cases the sequence of atoms in the h[] array is reversed (rev flag),
// as detailed in the comments.
//
// other parameters:
//
// h = hedron array data according to rev flag:
//   yes reversed                                     :  not reversed
//    0    1     2     3     4     5     6      7     :     0     1     2     3    4     5      6      7
//  len1  len3  atom1 atom3  a1    a2   a1-a2  a2-a3      len1  len3  atom1 atom3   a1    a3  a1-a2  a2-a3
//
// split: chop half of the hedron - to selectively print parts of a rotating
//   bond to be glued together.  top or bottom half selected by global caTop
//   (C-alpha top) variable, undef by default so bottom half.
//
// supporSel: enable support structure inside rotatable bond to print in place.
//  Please note the bond needs to be exactly parallel to the buildplate and the
//  layerHeight global variable above needs to be set correctly for the
//  structure to be correctly created by your slicer software.
//
 */

module hedron(h,rev=0,scal,split=0, supportSel) {

    newh = hFlip(h, rev);  // make a consistent hedron array regardless of rev flag

    bondRad = bondRadius * scal;
    difference() {
        union(){
            if (h[7]) {
                // central atom at 0,0,0
                atom(h[4],scal);
            }

            if (newh[5] && newh[7] != FemaleJoinBond) {  // not female join
                // comments for non-reversed case
                // atom 3 is len3 up on +z
                translate([0,0,newh[1]])
                    difference() {
                        atom(newh[3],scal * (newh[7] == SkinnyBond ? 0.7 : 1));  // if skinny bond make atom (C-beta) same diameter as bond
                        if (newh[7] == HBond) {  // make room for hbond magnet through atom - this branch not used for backbone N,O
                            translate([0,0,scal*hblen/2-magL-pClearance])
                                cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                        }
                    }
            }

            if (newh[7]) {
                // atom 2 - atom 3 bond from origin up +z distance len3
                bond(newh[1], bondRad, scal, newh[7], h[4], ver=1, supportSel=supportSel);
            }
            rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
                if (newh[6]) {
                    bond(newh[0], bondRad, scal, newh[6], h[4], ver=1, supportSel=supportSel);  // h[4] is center atom (atom 2)
                }
                if (newh[4] && newh[6] != FemaleJoinBond) {   // if draw atom 2 and atom1-atom2 not joiner
                    translate([0,0,newh[0]]) {
                        difference() {
                            atom(newh[2],scal * (newh[6] == SkinnyBond ? 0.7 : 1));  // put atom1 sphere len1 away on Z
                            if (newh[6] == HBond) {  // make room for hbond magnet through atom
                                translate([0,0,scal*hblen/2-magL-pClearance])
                                    cylinder(h=magL+pClearance,r=magR+pClearance,$fn=8);
                            }
                        }
                    }
                }
            }
        }

        if (split) {
            // top / bottom half cutter
            thick = 2*bondRadius * scal;
            Zdim = newh[0];
            Xdim = newh[1];

            cside = 7* defaultAtomRadius * atomScale * scal / 12 + (caTop ? pClearance : -pClearance);
            difference() {
                translate([-Xdim,((rev || caTop) ? 0 : -thick),-Zdim]) {
                    cube([2*Xdim,thick,2*Zdim]);
                }
                if (!caTop) {
                    rotate([0,(rev ? h[1] : 0),0])
                    rotate([45,0,0])
                    cube([cside, cside, cside],center=true);
                }
            }
            if (caTop) {
                //translate([tx+cside,0,tx+cside])
                    rotate([0,(rev ? h[1] : 0),0])
                        rotate([45,0,0])
                        cube([cside, cside, cside], center=true);
            }
        }

        if (newh[7] == FemaleJoinBond) {  // female join
            joiner(newh[1], scal, male=false, ver=1, supportSelect=supportSel);
        }

        if (newh[6] == FemaleJoinBond) {  // female join
            rotate([0, h[1], 0]) {                        // rotate following elements by angle2 about Y
            joiner(newh[0], scal, male=false, ver=1, supportSelect=supportSel);
            translate([0,0,newh[0]])
                atom(newh[2],scal+0.5,clearance);  // clearance for atom against join outer cylinder
            }
        }

        if (newh[7] == FemaleJoinBond || newh[6] == FemaleJoinBond) {  // female join both hedron arms
            translate([0,0,newh[1]]) atom(newh[3],scal+0.5,clearance);  // clearance for atom against join outer cylinder
        }
    }
}

/*
//
// Hook to call custom routines for specific hedra.
//
// Residue is h[h_residue]
// Sequence position is h[h_seqpos]
//
*/
module hedronDispatch(h,rev=0,scal) {
    // default action is just to pass to hedron()

    hedron(h, rev, scal, 0, (support ? 1 : 0));

    /*
    // Some examples for special handling for specific hedra below:
    // note use of h_seqpos, h_residue, h_class for selecting hedra

    // bool flag caTop (for rotatable bond part) needs to be a global variable
    // so hedron() above can see it.

caBase1 = false;   // only make bottom of N_C-alpha_C hedron
caBase2 = false;   // same as caBase1 but for case of reversed hedron (for testing, should be identical to caBase1 result)
amideOnly = false; // make only the first amide

    if (caTop) {
        // these examples select a specific sequence position (h[h_seqpos] == n)
        if (h[h_seqpos] == 1) {
            if (h[h_class] == "NCAC") {
                hedron(h, rev, scal, 1);
            } else if (h[h_class] == "CBCAC") {
                color("yellow") {  // ca-cb
                    hedron(h, rev, scal);
                }
            }
        }
    } else if (caBase1) {
        if (h[h_seqpos] == 1 && (h[h_class] == "NCAC")) {
            hedron(h, rev, scal, true, (support ? 1 : 0));
        }
    } else if (caBase2) {
        if (h[h_seqpos] == 5 && (h[h_class] == "NCAC")) {
            hedron(h, rev, scal, true, (support ? 1 : 0));
        }
    } else if (amideOnly) {
        if (h[h_seqpos] == 1) {
            if (h[h_class] == "CACN") {
                color("darkgray") {
                    hedron(h, rev, scal);
                }
            }  else if (h[h_class] == "CACO") {
                color("red") {   // c=o
                    hedron(h, rev, scal);
                }
            }  else if (h[h_class] == "CNCA") {
                color("cyan") {  // h=n
                    hedron(h, rev, scal);
                }
            }
        } else if ((h[h_seqpos] == 2) && (h[h_class] == "HNCA")) {
            color("cyan") {  // h=n
                hedron(h, rev, scal);
            }
        }
       // actions above select out only a single hedron
    } else {
        // actions below will process hedra all but handle selected ones differently

        if (h[h_class] == "NCAC") {
            if (h[h_seqpos] == 1) {
                if (! CCap && NCap) {  // make split rotatable bond for terminal NH3
                    hedron(h, rev, scal, true, (support ? 1 : 0));
                }
            } else if (h[h_seqpos] == 5) {  // make split rotatable bond for terminal COOH
                hedron(h, rev, scal, true, (support ? 2 : 0));  // note supportSel = 2
            } else {
                hedron(h, rev, scal, 0, (support ? 2 : 0));
            }
        } else if (h[h_class] == "CBCAC") {
            color("yellow") {                     // ca-cb -- color yellow in OpenSCAD renderer
                if (h[h_seqpos] == 1 ) {         // don't make here for N-term
                } else if (h[h_seqpos] == 5 ) {  // don't make here for C-term
                } else {
                    hedron(h, rev, scal);       // otherwise do make here
                }
            }
        } else if (h[h_class] == "HNCA") {
            color("cyan") { // color h-n in OenSCAD renderer
                if (h[h_seqpos] == 1) {
                    if (NCap) {                      // only make at N term if variable NCap is true
                        hedron(h, rev, scal, 0, (support ? 1 : 0));
                    }
                } else {
                    hedron(h, rev, scal, 0, (support ? 1 : 0));
                }
            }
        } else if (h[h_residue] == "P") {
            color("darkgray")   // highlight Prolines in OpenSCAD renderer
                hedron(h, rev, scal);
        } else {
            echo("unrecognised hedron", h[h_class]);
            color("pink")
                hedron(h, rev, scal, 0, (support ? 1 : 0));
        }
    }
    */
}

/*
//
// Generate a hedron rotated to specific angle d
//
*/
module d2(d,hedra,scal)
{
    tz = (d[d_reversed] ? hedra[d[d_h2ndx]][2] : hedra[d[d_h2ndx]][0]);      // get h2 len1 depending on reversed
    rotate(d[d_dangle1]) {                                                   // 4. rotate h2 to specified dihedral angle
        translate([0,0,tz]) {                                               // 3. translate h2 h2:len1 up +z
            rotate([180, 0, 0]) {                                          // 2. rotate h2r about X so h2:a3 in +z and h2:a1 in -z
                hedronDispatch(hedra[d[d_h2ndx]],(!d[d_reversed]),scal);  // 1. reverse hedron 2 orientation = h2r
            }
        }
    }
}

/*
//
// Generate two hedra at specified dihedral angle d
//
*/
module dihedron(d,hedra,scal)
{
    if (d[d_h1new])
        hedronDispatch(hedra[d[d_h1ndx]],d[d_reversed],scal);                // reverse h1 if dihedral reversed
    if (d[d_h2new])
        d2(d,hedra,scal);
}

/*
//
// Generate a residue consisting of the set of dihedra in the parameter 'r',
//   referring to hedra the table specified in the parameter 'hedra'.
//
*/
module residue(r,hedra, scal)
{
    for (d = r) {
        multmatrix(d[d_dihedralTransform]) {
            dihedron(d, hedra, scal);
        }
    }
}

/*
//
// Generate a chain of residues, each positioned by a supplied
// rotation/translation matrix.
//
*/
module chain(protein)
{
    chnD = protein[p_chainData];
    c = chnD[c_residues];
    dihedra = chnD[c_dihedra];
    hedra = chnD[c_hedra];
    for (r = c) {
        multmatrix(r[r_resTransform]) {
            residue(dihedra[r[r_resNdx]],hedra, protein[p_proteinScale]);
        }
    }
}

/*
//
// OpenSCAD array indices to reference protein data - tied to BioPython code
//
*/

// protein base level
p_pdbid = 0;
p_proteinScale = 1;
p_chainData = 2;

// chain level data
c_chainID = 0;
c_dihedra = 1;
c_hedra = 2;
c_residues = 3;

// hedra definitions
h_len1 = 0;
h_angle2 = 1;
h_len3 = 2;
h_atom1class = 3;
h_atom2class = 4;
h_atom3class = 5;
h_atom1state = 6;
h_atom2state = 7;
h_atom3state = 8;
h_bond1state = 9;
h_bond2state = 10;
h_residue = 11;
h_seqpos = 12;  // residue sequence position for first atom in hedra
h_class = 13;

// dihedra specifications for each residue in sequence, dihedral array
d_dangle1 = 0;
d_h1ndx = 1;
d_h2ndx = 2;
d_reversed = 3;
d_h1new = 4;
d_h2new = 5;
d_dihedralTransform = 6;

// residueSet: world transform for each residue in sequence array
r_resNdx = 0;
r_resID = 1;
r_resTransform = 2;


// use single default atom radius for all atoms if tubes = true, else use
// covalent radii from literature
atomData = ( tubes ?
            [   ["Csb","green" , defaultAtomRadius], ["Cres","green" , defaultAtomRadius], ["Cdb","green" , defaultAtomRadius],
                ["Osb","red" , defaultAtomRadius], ["Ores","red" , defaultAtomRadius], ["Odb","red" , defaultAtomRadius],
                ["Nsb","blue" , defaultAtomRadius], ["Nres","blue" , defaultAtomRadius], ["Ndb","blue" , defaultAtomRadius],
                ["Hsb","gray" , defaultAtomRadius],
                ["Ssb","yellow" , defaultAtomRadius] ]
            :

// covalent radii from Heyrovska, Raji : 'Atomic Structures of all the Twenty
// Essential Amino Acids and a Tripeptide, with Bond Lengths as Sums of Atomic
// Covalent Radii'  https://arxiv.org/pdf/0804.2488.pdf

            [   ["Csb","green" , 0.77], ["Cres","green" , 0.72], ["Cdb","green" , 0.67],
                ["Osb","red" , 0.67], ["Ores","red" , 0.635], ["Odb","red" , 0.60],
                ["Nsb","blue" , 0.70], ["Nres","blue" , 0.66], ["Ndb","blue" , 0.62],
                ["Hsb","gray" , 0.37],
                ["Ssb","yellow" , 1.04] ]
    );


// optionally include protein array data here [ write_SCAD(includeCode=False) ], e.g.:
// include <1rtm.scad>;
// or paste below

protein = [ "0PDB", 10.0,  // ID, protein_scale
 [
   "A", // chain id
   [  // residue array of dihedra
     [ // 0 : (' ', 260, ' ') L backbone
      [ 133.84135, 0, 1, 0, 1, 1,     // 260_L_N:260_L_CA:260_L_C:260_L_O [ 260_L_N:260_L_CA:260_L_C -- 260_L_CA:260_L_C:260_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.32066, 0, 2, 0, 0, 1,     // 260_L_N:260_L_CA:260_L_C:261_D_N [ 260_L_N:260_L_CA:260_L_C -- 260_L_CA:260_L_C:261_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.71709, 2, 3, 0, 0, 1,     // 260_L_CA:260_L_C:261_D_N:261_D_CA [ 260_L_CA:260_L_C:261_D_N -- 260_L_C:261_D_N:261_D_CA ] 
        [ [ 0.34235263327929577, 0.6860807516828992, 0.641940711169622, 0.0 ], [ -0.3228499694933943, 0.7275253962372915, -0.6053715347023565, 0.0 ], [ -0.882361927830551, 0.0, 0.4705713849302288, 15.276088445864564 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.44796, 3, 4, 0, 0, 1,     // 260_L_C:261_D_N:261_D_CA:261_D_C [ 260_L_C:261_D_N:261_D_CA -- 261_D_N:261_D_CA:261_D_C ] 
        [ [ -0.7092906910472809, -0.6991734171572244, 0.08979559195395095, 8.57595448422235 ], [ 0.6989726965294737, -0.7140876453216446, -0.038935899954993754, -8.08741155891488 ], [ 0.09134486904028485, 0.03514779565888809, 0.9951988481506268, 21.562648417090283 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 260, ' ') L sidechain

     ],
     [ // 1 : (' ', 261, ' ') D backbone
      [ 161.35560, 4, 5, 0, 0, 1,     // 261_D_N:261_D_CA:261_D_C:261_D_O [ 261_D_N:261_D_CA:261_D_C -- 261_D_CA:261_D_C:261_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -18.17688, 4, 6, 0, 0, 1,     // 261_D_N:261_D_CA:261_D_C:262_A_N [ 261_D_N:261_D_CA:261_D_C -- 261_D_CA:261_D_C:262_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.18298, 6, 7, 0, 0, 1,     // 261_D_CA:261_D_C:262_A_N:262_A_CA [ 261_D_CA:261_D_C:262_A_N -- 261_D_C:262_A_N:262_A_CA ] 
        [ [ 0.4409060700194232, 0.3119516338163659, 0.841598488341874, 0.0 ], [ -0.14476545713990605, 0.9500979834518648, -0.2763273136332601, 0.0 ], [ -0.885801783605735, 0.0, 0.4640637888921075, 15.451907467026146 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -85.28694, 7, 8, 0, 0, 1,     // 261_D_C:262_A_N:262_A_CA:262_A_C [ 261_D_C:262_A_N:262_A_CA -- 262_A_N:262_A_CA:262_A_C ] 
        [ [ -0.9224097110427717, -0.3550319457725995, 0.1520284264698962, 11.271672448928348 ], [ 0.36430407353422883, -0.9305333335181036, 0.0372861531678057, -3.700898957296082 ], [ 0.12822974296300318, 0.08977768482539435, 0.9876725673658423, 21.66719313621288 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 261, ' ') D sidechain

     ],
     [ // 2 : (' ', 262, ' ') A backbone
      [ 164.21670, 8, 9, 0, 0, 1,     // 262_A_N:262_A_CA:262_A_C:262_A_O [ 262_A_N:262_A_CA:262_A_C -- 262_A_CA:262_A_C:262_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -17.30560, 8, 10, 0, 0, 1,     // 262_A_N:262_A_CA:262_A_C:263_L_N [ 262_A_N:262_A_CA:262_A_C -- 262_A_CA:262_A_C:263_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.41221, 10, 11, 0, 0, 1,     // 262_A_CA:262_A_C:263_L_N:263_L_CA [ 262_A_CA:262_A_C:263_L_N -- 262_A_C:263_L_N:263_L_CA ] 
        [ [ 0.4728251308715136, 0.2974681387342013, 0.8294269721044439, 0.0 ], [ -0.14731929910770114, 0.9547317457998399, -0.258426619512011, 0.0 ], [ -0.8687539465963603, 0.0, 0.49524396036019286, 15.454488724614007 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -92.55599, 11, 12, 0, 0, 1,     // 262_A_C:263_L_N:263_L_CA:263_L_C [ 262_A_C:263_L_N:263_L_CA -- 263_L_N:263_L_CA:263_L_C ] 
        [ [ -0.9599635275583385, -0.2758166632052111, 0.048940719815856945, 11.104878957173051 ], [ 0.27315173724105735, -0.9604096093061889, -0.054786045618623165, -3.459974687959538 ], [ 0.062114041890247954, -0.039224362973496386, 0.9972979971650299, 22.08511977356237 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 262, ' ') A sidechain

     ],
     [ // 3 : (' ', 263, ' ') L backbone
      [ -63.84219, 12, 13, 0, 0, 1,     // 263_L_N:263_L_CA:263_L_C:263_L_O [ 263_L_N:263_L_CA:263_L_C -- 263_L_CA:263_L_C:263_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 119.99163, 12, 14, 0, 0, 1,     // 263_L_N:263_L_CA:263_L_C:264_S_N [ 263_L_N:263_L_CA:263_L_C -- 263_L_CA:263_L_C:264_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 167.47266, 14, 15, 0, 0, 1,     // 263_L_CA:263_L_C:264_S_N:264_S_CA [ 263_L_CA:263_L_C:264_S_N -- 263_L_C:264_S_N:264_S_CA ] 
        [ [ -0.22070238611477533, -0.8660984348254116, -0.4485130521581472, 0.0 ], [ 0.38239674005662094, -0.4998734851869743, 0.7771095366813464, 0.0 ], [ -0.8972531359417552, 0.0, 0.44151648898165324, 15.273497167853991 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.32005, 15, 16, 0, 0, 1,     // 263_L_C:264_S_N:264_S_CA:264_S_C [ 263_L_C:264_S_N:264_S_CA -- 264_S_N:264_S_CA:264_S_C ] 
        [ [ 0.39731936932324335, 0.8933504470179733, -0.20991021312305802, -5.971032496367331 ], [ -0.9142785840391806, 0.40502888888212785, -0.006802200983504502, 10.345621547542276 ], [ 0.078942951096923, 0.19461905863428808, 0.9776970044387026, 21.15138475820641 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 263, ' ') L sidechain

     ],
     [ // 4 : (' ', 264, ' ') S backbone
      [ -30.28810, 16, 17, 0, 0, 1,     // 264_S_N:264_S_CA:264_S_C:264_S_O [ 264_S_N:264_S_CA:264_S_C -- 264_S_CA:264_S_C:264_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 154.38066, 16, 18, 0, 0, 1,     // 264_S_N:264_S_CA:264_S_C:265_P_N [ 264_S_N:264_S_CA:264_S_C -- 264_S_CA:264_S_C:265_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.91934, 18, 19, 0, 0, 1,     // 264_S_CA:264_S_C:265_P_N:265_P_CA [ 264_S_CA:264_S_C:265_P_N -- 264_S_C:265_P_N:265_P_CA ] 
        [ [ -0.4295826292923246, -0.43239014562456096, -0.7927783590494031, 0.0 ], [ 0.20599983560952054, -0.9016866207096405, 0.38016484024969954, 0.0 ], [ -0.8792171701798958, 0.0, 0.47642120823999456, 15.413470044130525 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.55196, 19, 20, 0, 0, 1,     // 264_S_C:265_P_N:265_P_CA:265_P_C [ 264_S_C:265_P_N:265_P_CA -- 265_P_N:265_P_CA:265_P_C ] 
        [ [ 0.8668315441920852, 0.4820456898468638, -0.12741674494989386, -10.691104478887206 ], [ -0.49124741963925456, 0.8694171434816084, -0.052818588754854186, 5.126757031035107 ], [ 0.08531732937298987, 0.10837796602788534, 0.9904419063167306, 21.838303486051025 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 264, ' ') S sidechain

     ],
     [ // 5 : (' ', 265, ' ') P backbone
      [ 141.34270, 20, 21, 0, 0, 1,     // 265_P_N:265_P_CA:265_P_C:265_P_O [ 265_P_N:265_P_CA:265_P_C -- 265_P_CA:265_P_C:265_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.63052, 20, 22, 0, 0, 1,     // 265_P_N:265_P_CA:265_P_C:266_E_N [ 265_P_N:265_P_CA:265_P_C -- 265_P_CA:265_P_C:266_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.42129, 22, 23, 0, 0, 1,     // 265_P_CA:265_P_C:266_E_N:266_E_CA [ 265_P_CA:265_P_C:266_E_N -- 265_P_C:266_E_N:266_E_CA ] 
        [ [ 0.37623190642617627, 0.5825560595745892, 0.7204706725744308, 0.0 ], [ -0.26965887293853075, 0.8127905249527259, -0.5163871171830235, 0.0 ], [ -0.8864161803760389, 0.0, 0.46288913917649194, 15.382551126993805 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.47955, 23, 24, 0, 0, 1,     // 265_P_C:266_E_N:266_E_CA:266_E_C [ 265_P_C:266_E_N:266_E_CA -- 266_E_N:266_E_CA:266_E_C ] 
        [ [ -0.7864139509184924, -0.6049043344089899, 0.12507535334342956, 9.627178143007752 ], [ 0.6068658730481693, -0.7943735650798787, -0.02616201887788873, -6.9001431385281045 ], [ 0.11518207295517038, 0.05532978687380399, 0.991802250831409, 21.567835792908113 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 265, ' ') P sidechain

     ],
     [ // 6 : (' ', 266, ' ') E backbone
      [ 129.33020, 24, 25, 0, 0, 1,     // 266_E_N:266_E_CA:266_E_C:266_E_O [ 266_E_N:266_E_CA:266_E_C -- 266_E_CA:266_E_C:266_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.27428, 24, 26, 0, 0, 1,     // 266_E_N:266_E_CA:266_E_C:267_Q_N [ 266_E_N:266_E_CA:266_E_C -- 266_E_CA:266_E_C:267_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.86179, 26, 27, 0, 0, 1,     // 266_E_CA:266_E_C:267_Q_N:267_Q_CA [ 266_E_CA:266_E_C:267_Q_N -- 266_E_C:267_Q_N:267_Q_CA ] 
        [ [ 0.3067654821989282, 0.7463394493139631, 0.5906541842135627, 0.0 ], [ -0.3439949685218077, 0.6655654936951962, -0.6623367989429171, 0.0 ], [ -0.8874471255026631, 0.0, 0.4609095349817146, 15.343532069490612 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.04744, 27, 28, 0, 0, 1,     // 266_E_C:267_Q_N:267_Q_CA:267_Q_C [ 266_E_C:267_Q_N:267_Q_CA -- 267_Q_N:267_Q_CA:267_Q_C ] 
        [ [ -0.6471761921356176, -0.7572652397422321, 0.08781988960925452, 7.895659904688133 ], [ 0.7562400840368899, -0.6522676121196173, -0.05145772513092631, -8.853888191406115 ], [ 0.09624921624989764, 0.03311070609201888, 0.9948063980058474, 21.504810584581527 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 266, ' ') E sidechain

     ],
     [ // 7 : (' ', 267, ' ') Q backbone
      [ 135.73525, 28, 29, 0, 0, 1,     // 267_Q_N:267_Q_CA:267_Q_C:267_Q_O [ 267_Q_N:267_Q_CA:267_Q_C -- 267_Q_CA:267_Q_C:267_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.08278, 28, 30, 0, 0, 1,     // 267_Q_N:267_Q_CA:267_Q_C:268_L_N [ 267_Q_N:267_Q_CA:267_Q_C -- 267_Q_CA:267_Q_C:268_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.90563, 30, 31, 0, 0, 1,     // 267_Q_CA:267_Q_C:268_L_N:268_L_CA [ 267_Q_CA:267_Q_C:268_L_N -- 267_Q_C:268_L_N:268_L_CA ] 
        [ [ 0.366656782951396, 0.6438937120409814, 0.6715383020422726, 0.0 ], [ -0.3085653950978377, 0.7651149505754581, -0.5651429105580549, 0.0 ], [ -0.8776959612894708, 0.0, 0.47921790402295256, 15.297985904993082 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.19100, 31, 32, 0, 0, 1,     // 267_Q_C:268_L_N:268_L_CA:268_L_C [ 267_Q_C:268_L_N:268_L_CA -- 268_L_N:268_L_CA:268_L_C ] 
        [ [ -0.738620793628813, -0.6684294322586688, 0.08741520067754947, 8.99208629936776 ], [ 0.6713389533129313, -0.7411308587874311, 0.005390725155570826, -7.567422152629514 ], [ 0.061182783394027135, 0.06266693101914843, 0.9961573785164676, 21.714848109384913 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 267, ' ') Q sidechain

     ],
     [ // 8 : (' ', 268, ' ') L backbone
      [ 133.08127, 32, 33, 0, 0, 1,     // 268_L_N:268_L_CA:268_L_C:268_L_O [ 268_L_N:268_L_CA:268_L_C -- 268_L_CA:268_L_C:268_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.17416, 32, 34, 0, 0, 1,     // 268_L_N:268_L_CA:268_L_C:269_V_N [ 268_L_N:268_L_CA:268_L_C -- 268_L_CA:268_L_C:269_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.90526, 34, 35, 0, 0, 1,     // 268_L_CA:268_L_C:269_V_N:269_V_CA [ 268_L_CA:268_L_C:269_V_N -- 268_L_C:269_V_N:269_V_CA ] 
        [ [ 0.3529939696746897, 0.6842182612327701, 0.6381540788625495, 0.0 ], [ -0.33118392999217816, 0.7292772936240404, -0.598725173614181, 0.0 ], [ -0.8750499767946058, 0.0, 0.4840325795974479, 15.387386667077967 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.61549, 35, 36, 0, 0, 1,     // 268_L_C:269_V_N:269_V_CA:269_V_C [ 268_L_C:269_V_N:269_V_CA -- 269_V_N:269_V_CA:269_V_C ] 
        [ [ -0.7007482485947524, -0.7076775354344488, 0.09024631811221574, 8.54784051614913 ], [ 0.7104268446650693, -0.703767212809229, -0.0023256297991230196, -8.019704749329089 ], [ 0.06515819572870542, 0.06248372601048631, 0.9959167603335262, 21.870826062143554 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 268, ' ') L sidechain

     ],
     [ // 9 : (' ', 269, ' ') V backbone
      [ 142.14030, 36, 37, 0, 0, 1,     // 269_V_N:269_V_CA:269_V_C:269_V_O [ 269_V_N:269_V_CA:269_V_C -- 269_V_CA:269_V_C:269_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -36.94350, 36, 38, 0, 0, 1,     // 269_V_N:269_V_CA:269_V_C:270_L_N [ 269_V_N:269_V_CA:269_V_C -- 269_V_CA:269_V_C:270_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.74473, 38, 39, 0, 0, 1,     // 269_V_CA:269_V_C:270_L_N:270_L_CA [ 269_V_CA:269_V_C:270_L_N -- 269_V_C:270_L_N:270_L_CA ] 
        [ [ 0.38413850216995404, 0.601027215547228, 0.7008593991109557, 0.0 ], [ -0.28887568186742946, 0.7992285569044351, -0.5270527053854119, 0.0 ], [ -0.8769198661087864, 0.0, 0.4806366074528116, 15.322470515480383 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.37097, 39, 40, 0, 0, 1,     // 269_V_C:270_L_N:270_L_CA:270_L_C [ 269_V_C:270_L_N:270_L_CA -- 270_L_N:270_L_CA:270_L_C ] 
        [ [ -0.7666703932567673, -0.633685295581074, 0.10324463311882327, 9.369793121102623 ], [ 0.6384790870677508, -0.7694099378730342, 0.018783047658592827, -7.046170486752575 ], [ 0.06753490564622923, 0.08031994563332512, 0.9944786286556462, 21.748103936354994 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 269, ' ') V sidechain

     ],
     [ // 10 : (' ', 270, ' ') L backbone
      [ 129.56268, 40, 41, 0, 0, 1,     // 270_L_N:270_L_CA:270_L_C:270_L_O [ 270_L_N:270_L_CA:270_L_C -- 270_L_CA:270_L_C:270_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.37303, 40, 42, 0, 0, 1,     // 270_L_N:270_L_CA:270_L_C:271_T_N [ 270_L_N:270_L_CA:270_L_C -- 270_L_CA:270_L_C:271_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.25476, 42, 43, 0, 0, 1,     // 270_L_CA:270_L_C:271_T_N:271_T_CA [ 270_L_CA:270_L_C:271_T_N -- 270_L_C:271_T_N:271_T_CA ] 
        [ [ 0.3437850975182623, 0.6993270048624104, 0.626700524169662, 0.0 ], [ -0.33634242763182937, 0.7148018888266666, -0.6131329636421103, 0.0 ], [ -0.8767471574513866, 0.0, 0.4809515795804329, 15.358686762994965 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.09994, 43, 44, 0, 0, 1,     // 270_L_C:271_T_N:271_T_CA:271_T_C [ 270_L_C:271_T_N:271_T_CA -- 271_T_N:271_T_CA:271_T_C ] 
        [ [ -0.674615187616694, -0.7302289400210431, 0.10798168267181525, 8.375897419343996 ], [ 0.7353800595830342, -0.6775415957544104, 0.012391690369762577, -8.194566000545391 ], [ 0.06411331066592739, 0.08776719876076117, 0.9940756521599065, 21.786639092498024 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 270, ' ') L sidechain

     ],
     [ // 11 : (' ', 271, ' ') T backbone
      [ 138.42408, 44, 45, 0, 0, 1,     // 271_T_N:271_T_CA:271_T_C:271_T_O [ 271_T_N:271_T_CA:271_T_C -- 271_T_CA:271_T_C:271_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.53665, 44, 46, 0, 0, 1,     // 271_T_N:271_T_CA:271_T_C:272_L_N [ 271_T_N:271_T_CA:271_T_C -- 271_T_CA:271_T_C:272_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.91616, 46, 47, 0, 0, 1,     // 271_T_CA:271_T_C:272_L_N:272_L_CA [ 271_T_CA:271_T_C:272_L_N -- 271_T_C:272_L_N:272_L_CA ] 
        [ [ 0.3846001809231004, 0.6230150546831281, 0.6811276990932738, 0.0 ], [ -0.30632662612983713, 0.7822098450148648, -0.5425050750783186, 0.0 ], [ -0.8707736209588743, 0.0, 0.49168414764172597, 15.311266193517913 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.76270, 47, 48, 0, 0, 1,     // 271_T_C:272_L_N:272_L_CA:272_L_C [ 271_T_C:272_L_N:272_L_CA -- 272_L_N:272_L_CA:272_L_C ] 
        [ [ -0.7445911296680563, -0.6602674656219681, 0.09813726846000352, 9.130586803730882 ], [ 0.6657021744997771, -0.7453387498165539, 0.03620445950659766, -7.2723362829923754 ], [ 0.04924088226175613, 0.09228771241633169, 0.9945141093272857, 21.902342721513342 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 271, ' ') T sidechain

     ],
     [ // 12 : (' ', 272, ' ') L backbone
      [ 143.29759, 48, 49, 0, 0, 1,     // 272_L_N:272_L_CA:272_L_C:272_L_O [ 272_L_N:272_L_CA:272_L_C -- 272_L_CA:272_L_C:272_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.09677, 48, 50, 0, 0, 1,     // 272_L_N:272_L_CA:272_L_C:273_L_N [ 272_L_N:272_L_CA:272_L_C -- 272_L_CA:272_L_C:273_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.22809, 50, 51, 0, 0, 1,     // 272_L_CA:272_L_C:273_L_N:273_L_CA [ 272_L_CA:272_L_C:273_L_N -- 272_L_C:273_L_N:273_L_CA ] 
        [ [ 0.3902742201117288, 0.5749591256147889, 0.7191022437654184, 0.0 ], [ -0.2742564462222085, 0.8181821336794622, -0.5053329573986715, 0.0 ], [ -0.8789024034679543, 0.0, 0.47700164064524275, 15.401226195676976 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.46748, 51, 52, 0, 0, 1,     // 272_L_C:273_L_N:273_L_CA:273_L_C [ 272_L_C:273_L_N:273_L_CA -- 273_L_N:273_L_CA:273_L_C ] 
        [ [ -0.7931652566012403, -0.5993877516816013, 0.1078109403296519, 9.615866072833894 ], [ 0.6021034978928822, -0.7983679363513567, -0.00894516803821948, -6.757333998974207 ], [ 0.09143442210592594, 0.05781834778122865, 0.994131170979969, 21.779712610834054 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 272, ' ') L sidechain

     ],
     [ // 13 : (' ', 273, ' ') L backbone
      [ 132.52708, 52, 53, 0, 0, 1,     // 273_L_N:273_L_CA:273_L_C:273_L_O [ 273_L_N:273_L_CA:273_L_C -- 273_L_CA:273_L_C:273_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.74341, 52, 54, 0, 0, 1,     // 273_L_N:273_L_CA:273_L_C:274_E_N [ 273_L_N:273_L_CA:273_L_C -- 273_L_CA:273_L_C:274_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.27208, 54, 55, 0, 0, 1,     // 273_L_CA:273_L_C:274_E_N:274_E_CA [ 273_L_CA:273_L_C:274_E_N -- 273_L_C:274_E_N:274_E_CA ] 
        [ [ 0.3346882186236709, 0.7039330161296519, 0.6264679601680562, 0.0 ], [ -0.33170384930639696, 0.7102663646848352, -0.6208818305867131, 0.0 ], [ -0.8820183403250939, 0.0, 0.4712150754487456, 15.387409475822716 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.95570, 55, 56, 0, 0, 1,     // 273_L_C:274_E_N:274_E_CA:274_E_C [ 273_L_C:274_E_N:274_E_CA -- 274_E_N:274_E_CA:274_E_C ] 
        [ [ -0.6655431649079494, -0.7338217333618634, 0.13622760105830367, 8.364957986422295 ], [ 0.738839573868941, -0.6736144876496347, -0.018963283304366797, -8.290368793956473 ], [ 0.10568055511527061, 0.08802945912767082, 0.9904960548110315, 21.679342278962675 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 273, ' ') L sidechain

     ],
     [ // 14 : (' ', 274, ' ') E backbone
      [ 155.63670, 56, 57, 0, 0, 1,     // 274_E_N:274_E_CA:274_E_C:274_E_O [ 274_E_N:274_E_CA:274_E_C -- 274_E_CA:274_E_C:274_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -24.20305, 56, 58, 0, 0, 1,     // 274_E_N:274_E_CA:274_E_C:275_A_N [ 274_E_N:274_E_CA:274_E_C -- 274_E_CA:274_E_C:275_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.07531, 58, 59, 0, 0, 1,     // 274_E_CA:274_E_C:275_A_N:275_A_CA [ 274_E_CA:274_E_C:275_A_N -- 274_E_C:275_A_N:275_A_CA ] 
        [ [ 0.41548934101435725, 0.40997150888876754, 0.8119679608229148, 0.0 ], [ -0.18675485595484267, 0.9120983290739366, -0.3649647405953331, 0.0 ], [ -0.8902197657212185, 0.0, 0.45553130377533757, 15.408817281820127 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.04470, 59, 60, 0, 0, 1,     // 274_E_C:275_A_N:275_A_CA:275_A_C [ 274_E_C:275_A_N:275_A_CA -- 275_A_N:275_A_CA:275_A_C ] 
        [ [ -0.9040214817135129, -0.41662332776592065, 0.09576096992906222, 10.839576149547458 ], [ 0.4159183721236227, -0.9089656801648396, -0.028165582055165988, -4.87219113125295 ], [ 0.09877787368910214, 0.014366455503040056, 0.9950057972824776, 21.490050155554545 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 274, ' ') E sidechain

     ],
     [ // 15 : (' ', 275, ' ') A backbone
      [ 152.94334, 60, 61, 0, 0, 1,     // 275_A_N:275_A_CA:275_A_C:275_A_O [ 275_A_N:275_A_CA:275_A_C -- 275_A_CA:275_A_C:275_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -26.95629, 60, 62, 0, 0, 1,     // 275_A_N:275_A_CA:275_A_C:276_E_N [ 275_A_N:275_A_CA:275_A_C -- 275_A_CA:275_A_C:276_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -166.50966, 62, 63, 0, 0, 1,     // 275_A_CA:275_A_C:276_E_N:276_E_CA [ 275_A_CA:275_A_C:276_E_N -- 275_A_C:276_E_N:276_E_CA ] 
        [ [ 0.4320138160394385, 0.4533106591603446, 0.7796624327506464, 0.0 ], [ -0.2197070711409536, 0.8913525936976983, -0.39650896156209386, 0.0 ], [ -0.8746958703696424, 0.0, 0.4846721926810882, 15.326742261166004 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -50.56297, 63, 64, 0, 0, 1,     // 275_A_C:276_E_N:276_E_CA:276_E_C [ 275_A_C:276_E_N:276_E_CA -- 276_E_N:276_E_CA:276_E_C ] 
        [ [ -0.9403599444860433, -0.34002266571432144, -0.010380828807576698, 10.46074535007429 ], [ 0.334296997243791, -0.9180130984756402, -0.21330135644420825, -5.319968106312508 ], [ 0.06299755900025168, -0.20405033160394817, 0.9769313024631407, 21.829598077896776 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 275, ' ') A sidechain

     ],
     [ // 16 : (' ', 276, ' ') E backbone
      [ -38.97400, 64, 65, 0, 0, 1,     // 276_E_N:276_E_CA:276_E_C:276_E_O [ 276_E_N:276_E_CA:276_E_C -- 276_E_CA:276_E_C:276_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 139.38343, 64, 66, 0, 0, 1,     // 276_E_N:276_E_CA:276_E_C:277_P_N [ 276_E_N:276_E_CA:276_E_C -- 276_E_CA:276_E_C:277_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.77980, 66, 67, 0, 0, 1,     // 276_E_CA:276_E_C:277_P_N:277_P_CA [ 276_E_CA:276_E_C:277_P_N -- 276_E_C:277_P_N:277_P_CA ] 
        [ [ -0.35395850349450836, -0.650993748838262, -0.6715061554278069, 0.0 ], [ 0.3035567197699432, -0.7590830909548083, 0.5758872970547331, 0.0 ], [ -0.8846279984753141, 0.0, 0.4662974418904308, 15.398049080275559 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -73.30820, 67, 68, 0, 0, 1,     // 276_E_C:277_P_N:277_P_CA:277_P_C [ 276_E_C:277_P_N:277_P_CA -- 277_P_N:277_P_CA:277_P_C ] 
        [ [ 0.750361814673464, 0.6583836585437434, -0.05906018322441685, -9.036777150087394 ], [ -0.6581119594109879, 0.7524467313836427, 0.026693918976097804, 7.749988775209817 ], [ 0.062014481858485036, 0.018838115421341766, 0.9978974543745444, 21.67323514366329 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 276, ' ') E sidechain

     ],
     [ // 17 : (' ', 277, ' ') P backbone
      [ -31.30232, 68, 69, 0, 0, 1,     // 277_P_N:277_P_CA:277_P_C:277_P_O [ 277_P_N:277_P_CA:277_P_C -- 277_P_CA:277_P_C:277_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 149.86216, 68, 70, 0, 0, 1,     // 277_P_N:277_P_CA:277_P_C:278_P_N [ 277_P_N:277_P_CA:277_P_C -- 277_P_CA:277_P_C:278_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.70547, 70, 71, 0, 0, 1,     // 277_P_CA:277_P_C:278_P_N:278_P_CA [ 277_P_CA:277_P_C:278_P_N -- 277_P_C:278_P_N:278_P_CA ] 
        [ [ -0.4077315540745844, -0.5020820244833933, -0.762672026825805, 0.0 ], [ 0.23671363171630125, -0.8648200047932848, 0.4427787436960159, 0.0 ], [ -0.8818852739283062, 0.0, 0.47146406398409246, 15.403629428327072 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.61022, 71, 72, 0, 0, 1,     // 277_P_C:278_P_N:278_P_CA:278_P_C [ 277_P_C:278_P_N:278_P_CA -- 278_P_N:278_P_CA:278_P_C ] 
        [ [ 0.8460462725145036, 0.5246840291325343, -0.09440537239778325, -10.263391304417716 ], [ -0.5271249336492554, 0.8497871078498335, -0.0010842774015252079, 5.958539644812013 ], [ 0.07965556533966914, 0.05068077451524686, 0.9955332490705433, 21.748191709939327 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 277, ' ') P sidechain

     ],
     [ // 18 : (' ', 278, ' ') P backbone
      [ -22.70140, 72, 73, 0, 0, 1,     // 278_P_N:278_P_CA:278_P_C:278_P_O [ 278_P_N:278_P_CA:278_P_C -- 278_P_CA:278_P_C:278_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 160.75075, 72, 74, 0, 0, 1,     // 278_P_N:278_P_CA:278_P_C:279_H_N [ 278_P_N:278_P_CA:278_P_C -- 278_P_CA:278_P_C:279_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 170.66389, 74, 75, 0, 0, 1,     // 278_P_CA:278_P_C:279_H_N:279_H_CA [ 278_P_CA:278_P_C:279_H_N -- 278_P_C:279_H_N:279_H_CA ] 
        [ [ -0.4232140658809298, -0.32967835430720366, -0.8439206343856193, 0.0 ], [ 0.14778678644622795, -0.9440933124968601, 0.29469795223347434, 0.0 ], [ -0.8938953631116056, 0.0, 0.44827567389673395, 15.34915160577593 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.49134, 75, 76, 0, 0, 1,     // 278_P_C:279_H_N:279_H_CA:279_H_C [ 278_P_C:279_H_N:279_H_CA -- 279_H_N:279_H_CA:279_H_C ] 
        [ [ 0.9088231490153785, 0.39396756272920846, -0.1372226050291216, -11.25654368572566 ], [ -0.40822818428846075, 0.907612704809187, -0.09792306991459228, 3.9307966155194243 ], [ 0.08596646652223, 0.14501288765540374, 0.9856880992725455, 21.32842797474426 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 278, ' ') P sidechain

     ],
     [ // 19 : (' ', 279, ' ') H backbone
      [ -47.94328, 76, 77, 0, 0, 1,     // 279_H_N:279_H_CA:279_H_C:279_H_O [ 279_H_N:279_H_CA:279_H_C -- 279_H_CA:279_H_C:279_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 133.79802, 76, 78, 0, 0, 1,     // 279_H_N:279_H_CA:279_H_C:280_V_N [ 279_H_N:279_H_CA:279_H_C -- 279_H_CA:279_H_C:280_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.70177, 78, 79, 0, 0, 1,     // 279_H_CA:279_H_C:280_V_N:280_V_CA [ 279_H_CA:279_H_C:280_V_N -- 279_H_C:280_V_N:280_V_CA ] 
        [ [ -0.3056838659392952, -0.7217841164435644, -0.6209549608097111, 0.0 ], [ 0.31878621191041356, -0.6921182624738224, 0.647570584450858, 0.0 ], [ -0.897180430682823, 0.0, 0.44166421045833465, 15.366682598469675 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -97.73639, 79, 80, 0, 0, 1,     // 279_H_C:280_V_N:280_V_CA:280_V_C [ 279_H_C:280_V_N:280_V_CA -- 280_V_N:280_V_CA:280_V_C ] 
        [ [ 0.6818084346112027, 0.7233654492015866, -0.10899396953230212, -8.273331339312618 ], [ -0.7153473245647781, 0.6904495774011513, 0.10750621519049812, 8.627946226193488 ], [ 0.15302112184608904, 0.0046699002086186495, 0.9882118843147991, 21.251222676985073 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 279, ' ') H sidechain

     ],
     [ // 20 : (' ', 280, ' ') V backbone
      [ -40.40292, 80, 81, 0, 0, 1,     // 280_V_N:280_V_CA:280_V_C:280_V_O [ 280_V_N:280_V_CA:280_V_C -- 280_V_CA:280_V_C:280_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 141.52038, 80, 82, 0, 0, 1,     // 280_V_N:280_V_CA:280_V_C:281_L_N [ 280_V_N:280_V_CA:280_V_C -- 280_V_CA:280_V_C:281_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.84394, 82, 83, 0, 0, 1,     // 280_V_CA:280_V_C:281_L_N:281_L_CA [ 280_V_CA:280_V_C:281_L_N -- 280_V_C:281_L_N:281_L_CA ] 
        [ [ -0.3601745474135901, -0.6222361572428134, -0.6950514081815188, 0.0 ], [ 0.28628660645136456, -0.7828295884927299, 0.5524652155085104, 0.0 ], [ -0.887870640556369, 0.0, 0.4600931706057187, 15.296908467977255 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -113.26066, 83, 84, 0, 0, 1,     // 280_V_C:281_L_N:281_L_CA:281_L_C [ 280_V_C:281_L_N:281_L_CA -- 281_L_N:281_L_CA:281_L_C ] 
        [ [ 0.7283042004296711, 0.6622570474470645, -0.1760357768841155, -9.298365116007666 ], [ -0.671348956134966, 0.7410683736621008, 0.01040406912224208, 7.390853722794574 ], [ 0.1373447149801871, 0.11060410781026919, 0.9843287868403044, 21.452013309598282 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 280, ' ') V sidechain

     ],
     [ // 21 : (' ', 281, ' ') L backbone
      [ -32.47156, 84, 85, 0, 0, 1,     // 281_L_N:281_L_CA:281_L_C:281_L_O [ 281_L_N:281_L_CA:281_L_C -- 281_L_CA:281_L_C:281_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 148.85360, 84, 86, 0, 0, 1,     // 281_L_N:281_L_CA:281_L_C:282_I_N [ 281_L_N:281_L_CA:281_L_C -- 281_L_CA:281_L_C:282_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -172.28356, 86, 87, 0, 0, 1,     // 281_L_CA:281_L_C:282_I_N:282_I_CA [ 281_L_CA:281_L_C:282_I_N -- 281_L_C:282_I_N:282_I_CA ] 
        [ [ -0.40568939217762906, -0.5172265235435357, -0.7535866509019478, 0.0 ], [ 0.24517575734318994, -0.8558485399549783, 0.4554252130351645, 0.0 ], [ -0.8805140345762463, 0.0, 0.47402007859821815, 15.403039548267119 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -120.72675, 87, 88, 0, 0, 1,     // 281_L_C:282_I_N:282_I_CA:282_I_C [ 281_L_C:282_I_N:282_I_CA -- 282_I_N:282_I_CA:282_I_C ] 
        [ [ 0.8884798259280997, 0.45807076569354627, -0.02783473613494041, -10.066629958351516 ], [ -0.4496138702830279, 0.8810184631991815, 0.1471524214930517, 6.083702634383487 ], [ 0.09192913884014157, -0.11822707439110754, 0.9887221006496363, 21.735137812435866 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 281, ' ') L sidechain

     ],
     [ // 22 : (' ', 282, ' ') I backbone
      [ 150.80478, 88, 89, 0, 0, 1,     // 282_I_N:282_I_CA:282_I_C:282_I_O [ 282_I_N:282_I_CA:282_I_C -- 282_I_CA:282_I_C:282_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -26.88929, 88, 90, 0, 0, 1,     // 282_I_N:282_I_CA:282_I_C:283_S_N [ 282_I_N:282_I_CA:282_I_C -- 282_I_CA:282_I_C:283_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -152.61366, 90, 91, 0, 0, 1,     // 282_I_CA:282_I_C:283_S_N:283_S_CA [ 282_I_CA:282_I_C:283_S_N -- 282_I_C:283_S_N:283_S_CA ] 
        [ [ 0.3886172317317567, 0.4522680207014788, 0.8027641525703005, 0.0 ], [ -0.19706545398738526, 0.8918820759779662, -0.4070768593201768, 0.0 ], [ -0.9000788043532032, 0.0, 0.43572714622124264, 15.385976775301186 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  67.24087, 91, 92, 0, 0, 1,     // 282_I_C:283_S_N:283_S_CA:283_S_C [ 282_I_C:283_S_N:283_S_CA -- 283_S_N:283_S_CA:283_S_C ] 
        [ [ -0.9719138629751591, -0.22282085015030223, 0.07572523816403798, 10.735303956764453 ], [ 0.1681274047237809, -0.8825722418355174, -0.43908930039603006, -5.443807878783473 ], [ 0.16467124441610834, -0.4140254903743557, 0.8952464881711207, 21.212922761767555 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 282, ' ') I sidechain

     ],
     [ // 23 : (' ', 283, ' ') S backbone
      [ -71.89345, 92, 93, 0, 0, 1,     // 283_S_N:283_S_CA:283_S_C:283_S_O [ 283_S_N:283_S_CA:283_S_C -- 283_S_CA:283_S_C:283_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 108.95324, 92, 94, 0, 0, 1,     // 283_S_N:283_S_CA:283_S_C:284_R_N [ 283_S_N:283_S_CA:283_S_C -- 283_S_CA:283_S_C:284_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.17936, 94, 95, 0, 0, 1,     // 283_S_CA:283_S_C:284_R_N:284_R_CA [ 283_S_CA:283_S_C:284_R_N -- 283_S_C:284_R_N:284_R_CA ] 
        [ [ -0.1497980834083332, -0.9457839491895813, -0.288189617548151, 0.0 ], [ 0.4362012924877291, -0.32479643079221093, 0.8391875302795407, 0.0 ], [ -0.8872930556694472, 0.0, 0.4612060638811846, 15.342391985560024 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.62731, 95, 96, 0, 0, 1,     // 283_S_C:284_R_N:284_R_CA:284_R_C [ 283_S_C:284_R_N:284_R_CA -- 284_R_N:284_R_CA:284_R_C ] 
        [ [ 0.2867785999794867, 0.9536635124763941, -0.09101615002290418, -3.8714163612158328 ], [ -0.9507033700315507, 0.2950090641530917, 0.09556544500063428, 11.273287228363463 ], [ 0.11798786719174893, 0.059123236030373685, 0.9912534015865161, 21.538037466449865 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 283, ' ') S sidechain

     ],
     [ // 24 : (' ', 284, ' ') R backbone
      [ -49.28852, 96, 97, 0, 0, 1,     // 284_R_N:284_R_CA:284_R_C:284_R_O [ 284_R_N:284_R_CA:284_R_C -- 284_R_CA:284_R_C:284_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 133.42967, 96, 98, 0, 0, 1,     // 284_R_N:284_R_CA:284_R_C:285_P_N [ 284_R_N:284_R_CA:284_R_C -- 284_R_CA:284_R_C:285_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 168.53920, 98, 99, 0, 0, 1,     // 284_R_CA:284_R_C:285_P_N:285_P_CA [ 284_R_CA:284_R_C:285_P_N -- 284_R_C:285_P_N:285_P_CA ] 
        [ [ -0.34091037208049135, -0.7262187153206255, -0.5969811519018009, 0.0 ], [ 0.36012880725942326, -0.6874637281472094, 0.6306352865673193, 0.0 ], [ -0.8683820359668007, 0.0, 0.49589579511239446, 15.398100018599374 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.41082, 99, 100, 0, 0, 1,     // 284_R_C:285_P_N:285_P_CA:285_P_C [ 284_R_C:285_P_N:285_P_CA -- 285_P_N:285_P_CA:285_P_C ] 
        [ [ 0.6072277094711838, 0.7794765788110175, -0.1539180689050223, -8.042987443315946 ], [ -0.7940313879730452, 0.6021996136924275, -0.08288413709685682, 8.496401728956497 ], [ 0.028083158012927545, 0.17254532260760558, 0.9846011618326774, 22.079188090770607 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 284, ' ') R sidechain

     ],
     [ // 25 : (' ', 285, ' ') P backbone
      [ -27.40681, 100, 101, 0, 0, 1,     // 285_P_N:285_P_CA:285_P_C:285_P_O [ 285_P_N:285_P_CA:285_P_C -- 285_P_CA:285_P_C:285_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 156.09605, 100, 102, 0, 0, 1,     // 285_P_N:285_P_CA:285_P_C:286_S_N [ 285_P_N:285_P_CA:285_P_C -- 285_P_CA:285_P_C:286_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 166.30489, 102, 103, 0, 0, 1,     // 285_P_CA:285_P_C:286_S_N:286_S_CA [ 285_P_CA:285_P_C:286_S_N -- 285_P_C:286_S_N:286_S_CA ] 
        [ [ -0.40723969488790723, -0.4052046576322491, -0.8185139072372173, 0.0 ], [ 0.18049740492369024, -0.9142260034767835, 0.3627829948918246, 0.0 ], [ -0.8953080574435917, 0.0, 0.44544750788008947, 15.320933564529362 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.68136, 103, 104, 0, 0, 1,     // 285_P_C:286_S_N:286_S_CA:286_S_C [ 285_P_C:286_S_N:286_S_CA -- 286_S_N:286_S_CA:286_S_C ] 
        [ [ 0.8502391655749533, 0.490100542374574, -0.19208024283240435, -10.896017781955845 ], [ -0.5172894043649165, 0.8455002435146612, -0.1324424793955418, 4.829349786706692 ], [ 0.09749376110396477, 0.2119688575729864, 0.9724012391831913, 21.25070937392526 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 285, ' ') P sidechain

     ],
     [ // 26 : (' ', 286, ' ') S backbone
      [ -52.44855, 104, 105, 0, 0, 1,     // 286_S_N:286_S_CA:286_S_C:286_S_O [ 286_S_N:286_S_CA:286_S_C -- 286_S_CA:286_S_C:286_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 130.05734, 104, 106, 0, 0, 1,     // 286_S_N:286_S_CA:286_S_C:287_A_N [ 286_S_N:286_S_CA:286_S_C -- 286_S_CA:286_S_C:287_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 169.97194, 106, 107, 0, 0, 1,     // 286_S_CA:286_S_C:287_A_N:287_A_CA [ 286_S_CA:286_S_C:287_A_N -- 286_S_C:287_A_N:287_A_CA ] 
        [ [ -0.29549369300096556, -0.7654008284927067, -0.5717036375075196, 0.0 ], [ 0.3514408524537713, -0.6435538607938563, 0.6799468769572286, 0.0 ], [ -0.8883539861019467, 0.0, 0.45915922660530556, 15.299527658144418 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -70.61213, 107, 108, 0, 0, 1,     // 286_S_C:287_A_N:287_A_CA:287_A_C [ 286_S_C:287_A_N:287_A_CA -- 287_A_N:287_A_CA:287_A_C ] 
        [ [ 0.5664571388425947, 0.8051619422264317, -0.17561479620015422, -7.6153228856929225 ], [ -0.8196675321439504, 0.5725254629774915, -0.0189665756386062, 9.05716646429367 ], [ 0.08527277762162642, 0.15468949877922974, 0.984276237833731, 21.415713627834798 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 286, ' ') S sidechain

     ],
     [ // 27 : (' ', 287, ' ') A backbone
      [ -30.75100, 108, 109, 0, 0, 1,     // 287_A_N:287_A_CA:287_A_C:287_A_O [ 287_A_N:287_A_CA:287_A_C -- 287_A_CA:287_A_C:287_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 151.98110, 108, 110, 0, 0, 1,     // 287_A_N:287_A_CA:287_A_C:288_P_N [ 287_A_N:287_A_CA:287_A_C -- 287_A_CA:287_A_C:288_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -1.91747, 110, 111, 0, 0, 1,     // 287_A_CA:287_A_C:288_P_N:288_P_CA [ 287_A_CA:287_A_C:288_P_N -- 287_A_C:288_P_N:288_P_CA ] 
        [ [ -0.4353118696457881, -0.4697628432874578, -0.7680015932353081, 0.0 ], [ 0.23164368272590727, -0.8827926546287544, 0.4086787652749967, 0.0 ], [ -0.8699682640181006, 0.0, 0.4931077160229116, 15.498569969834529 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -75.18242, 111, 112, 0, 0, 1,     // 287_A_C:288_P_N:288_P_CA:288_P_C [ 287_A_C:288_P_N:288_P_CA -- 288_P_N:288_P_CA:288_P_C ] 
        [ [ 0.35801618190768353, -0.48406526815582473, -0.7984392460653934, -10.357448805026259 ], [ -0.16758958783058334, -0.8745475874430637, 0.45506070732188236, 5.511537249819496 ], [ -0.9185521996832948, -0.029108992815729186, -0.3942264869262772, 22.14873576253217 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 287, ' ') A sidechain

     ],
     [ // 28 : (' ', 288, ' ') P backbone
      [ -34.33996, 112, 113, 0, 0, 1,     // 288_P_N:288_P_CA:288_P_C:288_P_O [ 288_P_N:288_P_CA:288_P_C -- 288_P_CA:288_P_C:288_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 147.42468, 112, 114, 0, 0, 1,     // 288_P_N:288_P_CA:288_P_C:289_F_N [ 288_P_N:288_P_CA:288_P_C -- 288_P_CA:288_P_C:289_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.40958, 114, 115, 0, 0, 1,     // 288_P_CA:288_P_C:289_F_N:289_F_CA [ 288_P_CA:288_P_C:289_F_N -- 288_P_C:289_F_N:289_F_CA ] 
        [ [ -0.37935943813886935, -0.538407803626456, -0.7524649185770037, 0.0 ], [ 0.24238027459139008, -0.8426844231348622, 0.4807647714785961, 0.0 ], [ -0.8929379705130495, 0.0, 0.45017972057394623, 15.360618904985188 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.72011, 115, 116, 0, 0, 1,     // 288_P_C:289_F_N:289_F_CA:289_F_C [ 288_P_C:289_F_N:289_F_CA -- 289_F_N:289_F_CA:289_F_C ] 
        [ [ 0.815698333191398, 0.567041756701325, -0.11445468704186254, -10.041037585206633 ], [ -0.5714790755092255, 0.8205830584582237, -0.007423639704928991, 6.415418208717261 ], [ 0.08971006344828975, 0.0714639092719054, 0.9934007319242727, 21.367904271762363 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 288, ' ') P sidechain

     ],
     [ // 29 : (' ', 289, ' ') F backbone
      [ -48.40168, 116, 117, 0, 0, 1,     // 289_F_N:289_F_CA:289_F_C:289_F_O [ 289_F_N:289_F_CA:289_F_C -- 289_F_CA:289_F_C:289_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 134.92181, 116, 118, 0, 0, 1,     // 289_F_N:289_F_CA:289_F_C:290_T_N [ 289_F_N:289_F_CA:289_F_C -- 289_F_CA:289_F_C:290_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 160.72887, 118, 119, 0, 0, 1,     // 289_F_CA:289_F_C:290_T_N:290_T_CA [ 289_F_CA:289_F_C:290_T_N -- 289_F_C:290_T_N:290_T_CA ] 
        [ [ -0.3096424154498603, -0.7080711106322335, -0.6346312920447679, 0.0 ], [ 0.31048870855204286, -0.7061411348222361, 0.6363658221287888, 0.0 ], [ -0.8987315152013201, 0.0, 0.43849933133807545, 15.337224820816616 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -94.81319, 119, 120, 0, 0, 1,     // 289_F_C:290_T_N:290_T_CA:290_T_C [ 289_F_C:290_T_N:290_T_CA -- 290_T_N:290_T_CA:290_T_C ] 
        [ [ 0.5554936834118359, 0.7705900304889836, -0.31243843009549666, -8.442297217320812 ], [ -0.822623212830526, 0.5641008931101, -0.0712827616241877, 8.465371116583654 ], [ 0.12131701200546546, 0.296616228995271, 0.9472597295856563, 21.170441258702773 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 289, ' ') F sidechain

     ],
     [ // 30 : (' ', 290, ' ') T backbone
      [ -26.76335, 120, 121, 0, 0, 1,     // 290_T_N:290_T_CA:290_T_C:290_T_O [ 290_T_N:290_T_CA:290_T_C -- 290_T_CA:290_T_C:290_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 157.93286, 120, 122, 0, 0, 1,     // 290_T_N:290_T_CA:290_T_C:291_E_N [ 290_T_N:290_T_CA:290_T_C -- 290_T_CA:290_T_C:291_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.33892, 122, 123, 0, 0, 1,     // 290_T_CA:290_T_C:291_E_N:291_E_CA [ 290_T_CA:290_T_C:291_E_N -- 290_T_C:291_E_N:291_E_CA ] 
        [ [ -0.4313759939538774, -0.37569275314720957, -0.8202254001632564, 0.0 ], [ 0.17487546328132228, -0.9267442771512915, 0.3325110781755595, 0.0 ], [ -0.8850611979871489, 0.0, 0.46547467795525976, 15.29603813636479 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.77235, 123, 124, 0, 0, 1,     // 290_T_C:291_E_N:291_E_CA:291_E_C [ 290_T_C:291_E_N:291_E_CA -- 291_E_N:291_E_CA:291_E_C ] 
        [ [ 0.9062375963167553, 0.4024713222322539, -0.12942277157691853, -10.97946623424583 ], [ -0.4063968853723933, 0.9136864228362479, -0.00432345744966673, 4.4509645209896975 ], [ 0.11651176155928525, 0.05651509095208615, 0.9915800794252725, 21.52684189895128 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 290, ' ') T sidechain

     ],
     [ // 31 : (' ', 291, ' ') E backbone
      [ 122.79676, 124, 125, 0, 0, 1,     // 291_E_N:291_E_CA:291_E_C:291_E_O [ 291_E_N:291_E_CA:291_E_C -- 291_E_CA:291_E_C:291_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -51.90846, 124, 126, 0, 0, 1,     // 291_E_N:291_E_CA:291_E_C:292_A_N [ 291_E_N:291_E_CA:291_E_C -- 291_E_CA:291_E_C:292_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.07283, 126, 127, 0, 0, 1,     // 291_E_CA:291_E_C:292_A_N:292_A_CA [ 291_E_CA:291_E_C:292_A_N -- 291_E_C:292_A_N:292_A_CA ] 
        [ [ 0.27690216208219137, 0.7870261510088417, 0.551284890290325, 0.0 ], [ -0.35325386000610165, 0.6169196362802921, -0.7032935893086055, 0.0 ], [ -0.893608920627473, 0.0, 0.4488464068865904, 15.380773510875379 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.40480, 127, 128, 0, 0, 1,     // 291_E_C:292_A_N:292_A_CA:292_A_C [ 291_E_C:292_A_N:292_A_CA -- 292_A_N:292_A_CA:292_A_C ] 
        [ [ -0.5971005657100834, -0.7958929848405727, 0.10012627582339328, 7.373964276979281 ], [ 0.7921748464407645, -0.6046910566758517, -0.08250902158426644, -9.40722645429125 ], [ 0.12621381499364254, 0.030051333710931244, 0.9915477750703416, 21.384524799050872 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 291, ' ') E sidechain

     ],
     [ // 32 : (' ', 292, ' ') A backbone
      [ 135.75068, 128, 129, 0, 0, 1,     // 292_A_N:292_A_CA:292_A_C:292_A_O [ 292_A_N:292_A_CA:292_A_C -- 292_A_CA:292_A_C:292_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.61660, 128, 130, 0, 0, 1,     // 292_A_N:292_A_CA:292_A_C:293_S_N [ 292_A_N:292_A_CA:292_A_C -- 292_A_CA:292_A_C:293_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.16698, 130, 131, 0, 0, 1,     // 292_A_CA:292_A_C:293_S_N:293_S_CA [ 292_A_CA:292_A_C:293_S_N -- 292_A_C:293_S_N:293_S_CA ] 
        [ [ 0.3621011213335059, 0.6770891955265449, 0.6406504501131908, 0.0 ], [ -0.3331627088941655, 0.7359009588940731, -0.5894509208586286, 0.0 ], [ -0.8705661303607666, 0.0, 0.4920514329505411, 15.401946644007813 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.50199, 131, 132, 0, 0, 1,     // 292_A_C:293_S_N:293_S_CA:293_S_C [ 292_A_C:293_S_N:293_S_CA -- 293_S_N:293_S_CA:293_S_C ] 
        [ [ -0.7297721523836116, -0.6822820298421228, 0.04386157042222112, 8.59153455318201 ], [ 0.6819895974367299, -0.7309795070814417, -0.02364633618680701, -7.904915938273514 ], [ 0.04819537937891706, 0.012656697100160728, 0.9987577351014795, 22.00067280956405 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 292, ' ') A sidechain

     ],
     [ // 33 : (' ', 293, ' ') S backbone
      [ 140.23190, 132, 133, 0, 0, 1,     // 293_S_N:293_S_CA:293_S_C:293_S_O [ 293_S_N:293_S_CA:293_S_C -- 293_S_CA:293_S_C:293_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.13576, 132, 134, 0, 0, 1,     // 293_S_N:293_S_CA:293_S_C:294_M_N [ 293_S_N:293_S_CA:293_S_C -- 293_S_CA:293_S_C:294_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.65495, 134, 135, 0, 0, 1,     // 293_S_CA:293_S_C:294_M_N:294_M_CA [ 293_S_CA:293_S_C:294_M_N -- 293_S_C:294_M_N:294_M_CA ] 
        [ [ 0.3584342563017519, 0.6037056043261626, 0.712084564651273, 0.0 ], [ -0.27143348837124376, 0.797207340222845, -0.5392442100620414, 0.0 ], [ -0.8932237935142726, 0.0, 0.4496123382425932, 15.331670271628024 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.53777, 135, 136, 0, 0, 1,     // 293_S_C:294_M_N:294_M_CA:294_M_C [ 293_S_C:294_M_N:294_M_CA -- 294_M_N:294_M_CA:294_M_C ] 
        [ [ -0.7694537626921297, -0.623591324870998, 0.13811142829091694, 9.519424500026487 ], [ 0.6241952193806737, -0.7800112161442432, -0.04430384623812205, -7.2088271528198895 ], [ 0.13535595731712807, 0.05211873209137051, 0.9894252890359897, 21.34226356440189 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 293, ' ') S sidechain

     ],
     [ // 34 : (' ', 294, ' ') M backbone
      [ 117.75519, 136, 137, 0, 0, 1,     // 294_M_N:294_M_CA:294_M_C:294_M_O [ 294_M_N:294_M_CA:294_M_C -- 294_M_CA:294_M_C:294_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.66253, 136, 138, 0, 0, 1,     // 294_M_N:294_M_CA:294_M_C:295_M_N [ 294_M_N:294_M_CA:294_M_C -- 294_M_CA:294_M_C:295_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.06424, 138, 139, 0, 0, 1,     // 294_M_CA:294_M_C:295_M_N:295_M_CA [ 294_M_CA:294_M_C:295_M_N -- 294_M_C:295_M_N:295_M_CA ] 
        [ [ 0.24352485045975505, 0.8449122387469163, 0.47625513753053383, 0.0 ], [ -0.38466109803747556, 0.5349049530670601, -0.7522715805086151, 0.0 ], [ -0.8903546972219316, 0.0, 0.4552675182075287, 15.254738242556373 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.09368, 139, 140, 0, 0, 1,     // 294_M_C:295_M_N:295_M_CA:295_M_C [ 294_M_C:295_M_N:295_M_CA -- 295_M_N:295_M_CA:295_M_C ] 
        [ [ -0.5555252925239134, -0.8313309324405984, 0.01674903381289757, 6.355337672837356 ], [ 0.8256790345584941, -0.5539038388494595, -0.1069797606951198, -10.038610692162422 ], [ 0.09821293833692939, -0.04560063678587384, 0.9941201138031306, 21.33000878910869 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 294, ' ') M sidechain

     ],
     [ // 35 : (' ', 295, ' ') M backbone
      [ 137.88080, 140, 141, 0, 0, 1,     // 295_M_N:295_M_CA:295_M_C:295_M_O [ 295_M_N:295_M_CA:295_M_C -- 295_M_CA:295_M_C:295_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -36.37946, 140, 142, 0, 0, 1,     // 295_M_N:295_M_CA:295_M_C:296_M_N [ 295_M_N:295_M_CA:295_M_C -- 295_M_CA:295_M_C:296_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.37344, 142, 143, 0, 0, 1,     // 295_M_CA:295_M_C:296_M_N:296_M_CA [ 295_M_CA:295_M_C:296_M_N -- 295_M_C:296_M_N:296_M_CA ] 
        [ [ 0.39847239918510885, 0.5931302773107471, 0.6995828908892352, 0.0 ], [ -0.2935587341939634, 0.8051064986306322, -0.5153899450330417, 0.0 ], [ -0.8689321128063465, 0.0, 0.494931291528328, 15.335739024004006 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.37930, 143, 144, 0, 0, 1,     // 295_M_C:296_M_N:296_M_CA:296_M_C [ 295_M_C:296_M_N:296_M_CA -- 296_M_N:296_M_CA:296_M_C ] 
        [ [ -0.7699828692354416, -0.6293405320849859, 0.10515167977233406, 9.36670121031274 ], [ 0.6346279286154799, -0.7724457997187617, 0.023976628580094573, -6.900545574789287 ], [ 0.06613450918531091, 0.08519378599304835, 0.9941671114671762, 21.962364099822267 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 295, ' ') M sidechain

     ],
     [ // 36 : (' ', 296, ' ') M backbone
      [ 145.20214, 144, 145, 0, 0, 1,     // 296_M_N:296_M_CA:296_M_C:296_M_O [ 296_M_N:296_M_CA:296_M_C -- 296_M_CA:296_M_C:296_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.28284, 144, 146, 0, 0, 1,     // 296_M_N:296_M_CA:296_M_C:297_S_N [ 296_M_N:296_M_CA:296_M_C -- 296_M_CA:296_M_C:297_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.74503, 146, 147, 0, 0, 1,     // 296_M_CA:296_M_C:297_S_N:297_S_CA [ 296_M_CA:296_M_C:297_S_N -- 296_M_C:297_S_N:297_S_CA ] 
        [ [ 0.3949600856530694, 0.5487724960298348, 0.7367871323130617, 0.0 ], [ -0.2592710040894125, 0.8359717385182258, -0.48366289793334666, 0.0 ], [ -0.8813541156535142, 0.0, 0.4724562655110122, 15.330031634484707 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.24698, 147, 148, 0, 0, 1,     // 296_M_C:297_S_N:297_S_CA:297_S_C [ 296_M_C:297_S_N:297_S_CA -- 297_S_N:297_S_CA:297_S_C ] 
        [ [ -0.8058786578915909, -0.5826396097596391, 0.10531226848750218, 9.823057168579787 ], [ 0.5877407521552468, -0.8087120178171109, 0.02335980509682563, -6.4483323450626004 ], [ 0.07155694942374952, 0.08072148027205882, 0.9941647980148234, 21.62895394492 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 296, ' ') M sidechain

     ],
     [ // 37 : (' ', 297, ' ') S backbone
      [ 137.24576, 148, 149, 0, 0, 1,     // 297_S_N:297_S_CA:297_S_C:297_S_O [ 297_S_N:297_S_CA:297_S_C -- 297_S_CA:297_S_C:297_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.06344, 148, 150, 0, 0, 1,     // 297_S_N:297_S_CA:297_S_C:298_L_N [ 297_S_N:297_S_CA:297_S_C -- 297_S_CA:297_S_C:298_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.60640, 150, 151, 0, 0, 1,     // 297_S_CA:297_S_C:298_L_N:298_L_CA [ 297_S_CA:297_S_C:298_L_N -- 297_S_C:298_L_N:298_L_CA ] 
        [ [ 0.3503349316818685, 0.6436354256462361, 0.6804402064080642, 0.0 ], [ -0.2946275628312674, 0.7653322408295553, -0.5722422217618333, 0.0 ], [ -0.8890781938972343, 0.0, 0.45775535511496956, 15.268071753037848 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.74102, 151, 152, 0, 0, 1,     // 297_S_C:298_L_N:298_L_CA:298_L_C [ 297_S_C:298_L_N:298_L_CA -- 298_L_N:298_L_CA:298_L_C ] 
        [ [ -0.751129236830366, -0.651965385345506, 0.10366294366801425, 9.088382083198235 ], [ 0.6494798148649267, -0.7579403588830917, -0.06084720584538586, -7.6432225881591584 ], [ 0.11824060073279677, 0.021622874171950347, 0.9927495211033301, 21.382136631532113 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 297, ' ') S sidechain

     ],
     [ // 38 : (' ', 298, ' ') L backbone
      [ 133.73919, 152, 153, 0, 0, 1,     // 298_L_N:298_L_CA:298_L_C:298_L_O [ 298_L_N:298_L_CA:298_L_C -- 298_L_CA:298_L_C:298_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.29273, 152, 154, 0, 0, 1,     // 298_L_N:298_L_CA:298_L_C:299_T_N [ 298_L_N:298_L_CA:298_L_C -- 298_L_CA:298_L_C:299_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.43280, 154, 155, 0, 0, 1,     // 298_L_CA:298_L_C:299_T_N:299_T_CA [ 298_L_CA:298_L_C:299_T_N -- 298_L_C:299_T_N:299_T_CA ] 
        [ [ 0.3695809868061474, 0.6729186637905109, 0.6407732563971342, 0.0 ], [ -0.33620711484645255, 0.7397164807697563, -0.58291020235044, 0.0 ], [ -0.8662416926689794, 0.0, 0.49962518939899486, 15.382091819301483 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.53061, 155, 156, 0, 0, 1,     // 298_L_C:299_T_N:299_T_CA:299_T_C [ 298_L_C:299_T_N:299_T_CA -- 299_T_N:299_T_CA:299_T_C ] 
        [ [ -0.7097553978960878, -0.7002109877889223, 0.0771482192732154, 8.581997829256013 ], [ 0.7032945361488091, -0.7105960189667988, 0.020738689731235275, -7.807026965123871 ], [ 0.040299759063768294, 0.06897731807049905, 0.9968039220484656, 22.073666800087423 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 298, ' ') L sidechain

     ],
     [ // 39 : (' ', 299, ' ') T backbone
      [ 136.55737, 156, 157, 0, 0, 1,     // 299_T_N:299_T_CA:299_T_C:299_T_O [ 299_T_N:299_T_CA:299_T_C -- 299_T_CA:299_T_C:299_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.43586, 156, 158, 0, 0, 1,     // 299_T_N:299_T_CA:299_T_C:300_K_N [ 299_T_N:299_T_CA:299_T_C -- 299_T_CA:299_T_C:300_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.83154, 158, 159, 0, 0, 1,     // 299_T_CA:299_T_C:300_K_N:300_K_CA [ 299_T_CA:299_T_C:300_K_N -- 299_T_C:300_K_N:300_K_CA ] 
        [ [ 0.36607366161324206, 0.6352140516650814, 0.6800685133428137, 0.0 ], [ -0.3010802196754518, 0.7723361370331129, -0.5593278043803565, 0.0 ], [ -0.880534369342421, 0.0, 0.47398230389619483, 15.276600794347207 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.49286, 159, 160, 0, 0, 1,     // 299_T_C:300_K_N:300_K_CA:300_K_C [ 299_T_C:300_K_N:300_K_CA -- 300_K_N:300_K_CA:300_K_C ] 
        [ [ -0.7354155514208187, -0.6656088636726064, 0.12699924144999755, 9.110231976972765 ], [ 0.6702723017573333, -0.7420733148572192, -0.007889035032557183, -7.4927833727059525 ], [ 0.0994937597306458, 0.07932235483947898, 0.9918714411643176, 21.626091554107827 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 299, ' ') T sidechain

     ],
     [ // 40 : (' ', 300, ' ') K backbone
      [ 147.78576, 160, 161, 0, 0, 1,     // 300_K_N:300_K_CA:300_K_C:300_K_O [ 300_K_N:300_K_CA:300_K_C -- 300_K_CA:300_K_C:300_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -32.31945, 160, 162, 0, 0, 1,     // 300_K_N:300_K_CA:300_K_C:301_L_N [ 300_K_N:300_K_CA:300_K_C -- 300_K_CA:300_K_C:301_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 168.79262, 162, 163, 0, 0, 1,     // 300_K_CA:300_K_C:301_L_N:301_L_CA [ 300_K_CA:300_K_C:301_L_N -- 300_K_C:301_L_N:301_L_CA ] 
        [ [ 0.4072088995628496, 0.5346392553462465, 0.7405010322475074, 0.0 ], [ -0.2576202989019654, 0.8450803906391456, -0.468477229917155, 0.0 ], [ -0.8762492189499943, 0.0, 0.48185818068133385, 15.308004450843281 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.00087, 163, 164, 0, 0, 1,     // 300_K_C:301_L_N:301_L_CA:301_L_C [ 300_K_C:301_L_N:301_L_CA -- 301_L_N:301_L_CA:301_L_C ] 
        [ [ -0.7837177746715038, -0.6035891382879253, 0.14651485182323987, 9.89641075489931 ], [ 0.6188698474031414, -0.7788935865445489, 0.10161147974014434, -6.260954265662574 ], [ 0.05278789292213177, 0.17030834677316686, 0.9839758662590421, 21.747788156691072 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 300, ' ') K sidechain

     ],
     [ // 41 : (' ', 301, ' ') L backbone
      [ 126.84112, 164, 165, 0, 0, 1,     // 301_L_N:301_L_CA:301_L_C:301_L_O [ 301_L_N:301_L_CA:301_L_C -- 301_L_CA:301_L_C:301_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -52.57413, 164, 166, 0, 0, 1,     // 301_L_N:301_L_CA:301_L_C:302_A_N [ 301_L_N:301_L_CA:301_L_C -- 301_L_CA:301_L_C:302_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -175.21833, 166, 167, 0, 0, 1,     // 301_L_CA:301_L_C:302_A_N:302_A_CA [ 301_L_CA:301_L_C:302_A_N -- 301_L_C:302_A_N:302_A_CA ] 
        [ [ 0.28034369728496183, 0.7941402536127118, 0.5392111543584159, 0.0 ], [ -0.3663313573300086, 0.6077345288791297, -0.7045992329294829, 0.0 ], [ -0.8872478503943254, 0.0, 0.46129302181005166, 15.276432417294167 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.30984, 167, 168, 0, 0, 1,     // 301_L_C:302_A_N:302_A_CA:302_A_C [ 301_L_C:302_A_N:302_A_CA -- 302_A_N:302_A_CA:302_A_C ] 
        [ [ -0.6403979153056486, -0.768007126777986, 0.007454078776216083, 7.203281294341249 ], [ 0.7614729568759976, -0.6361563892218152, -0.12435426972479662, -9.412688208587168 ], [ 0.10024692523316832, -0.0739601356846042, 0.9922098831954913, 21.438811233749064 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 301, ' ') L sidechain

     ],
     [ // 42 : (' ', 302, ' ') A backbone
      [ 133.76884, 168, 169, 0, 0, 1,     // 302_A_N:302_A_CA:302_A_C:302_A_O [ 302_A_N:302_A_CA:302_A_C -- 302_A_CA:302_A_C:302_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.78276, 168, 170, 0, 0, 1,     // 302_A_N:302_A_CA:302_A_C:303_D_N [ 302_A_N:302_A_CA:302_A_C -- 302_A_CA:302_A_C:303_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.95256, 170, 171, 0, 0, 1,     // 302_A_CA:302_A_C:303_D_N:303_D_CA [ 302_A_CA:302_A_C:303_D_N -- 302_A_C:303_D_N:303_D_CA ] 
        [ [ 0.35624972096005764, 0.6663080700315894, 0.6550722800780533, 0.0 ], [ -0.31833112598876706, 0.7456765758764174, -0.585347536439631, 0.0 ], [ -0.8784938420629963, 0.0, 0.4777536702709831, 15.356036258628338 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.34237, 171, 172, 0, 0, 1,     // 302_A_C:303_D_N:303_D_CA:303_D_C [ 302_A_C:303_D_N:303_D_CA -- 303_D_N:303_D_CA:303_D_C ] 
        [ [ -0.709141982801798, -0.6950672921675711, 0.11831782446772225, 8.74860355108374 ], [ 0.6993091337631557, -0.7147780511446202, -0.007685898602208874, -7.817417545593134 ], [ 0.08991320071802161, 0.07729034196288, 0.9929460304447074, 21.736519120424894 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 302, ' ') A sidechain

     ],
     [ // 43 : (' ', 303, ' ') D backbone
      [ 126.97594, 172, 173, 0, 0, 1,     // 303_D_N:303_D_CA:303_D_C:303_D_O [ 303_D_N:303_D_CA:303_D_C -- 303_D_CA:303_D_C:303_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.57881, 172, 174, 0, 0, 1,     // 303_D_N:303_D_CA:303_D_C:304_K_N [ 303_D_N:303_D_CA:303_D_C -- 303_D_CA:303_D_C:304_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.53092, 174, 175, 0, 0, 1,     // 303_D_CA:303_D_C:304_K_N:304_K_CA [ 303_D_CA:303_D_C:304_K_N -- 303_D_C:304_K_N:304_K_CA ] 
        [ [ 0.32197286121988156, 0.7263204883998726, 0.60728249173548, 0.0 ], [ -0.34022459564030455, 0.6873562017837407, -0.6417076253176526, 0.0 ], [ -0.8835047827597053, 0.0, 0.46842213743665645, 15.450262384112728 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.85818, 175, 176, 0, 0, 1,     // 303_D_C:304_K_N:304_K_CA:304_K_C [ 303_D_C:304_K_N:304_K_CA -- 304_K_N:304_K_CA:304_K_C ] 
        [ [ -0.6468363656927719, -0.7537008974316305, 0.1163515072017617, 8.142663045167025 ], [ 0.7583645471316429, -0.6518008688629047, -0.006232255079889595, -8.604247673176673 ], [ 0.0805352697343524, 0.084205608841202, 0.9931885449241235, 21.731035632067652 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 303, ' ') D sidechain

     ],
     [ // 44 : (' ', 304, ' ') K backbone
      [ 135.54739, 176, 177, 0, 0, 1,     // 304_K_N:304_K_CA:304_K_C:304_K_O [ 304_K_N:304_K_CA:304_K_C -- 304_K_CA:304_K_C:304_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.39194, 176, 178, 0, 0, 1,     // 304_K_N:304_K_CA:304_K_C:305_E_N [ 304_K_N:304_K_CA:304_K_C -- 304_K_CA:304_K_C:305_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.39681, 178, 179, 0, 0, 1,     // 304_K_CA:304_K_C:305_E_N:305_E_CA [ 304_K_CA:304_K_C:305_E_N -- 304_K_C:305_E_N:305_E_CA ] 
        [ [ 0.3329056267688649, 0.7119272565392512, 0.6183311613222517, 0.0 ], [ -0.3374916392950391, 0.7022532174340643, -0.6268491142288166, 0.0 ], [ -0.8804960176352739, 0.0, 0.47405354436859093, 15.314645231785365 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.93528, 179, 180, 0, 0, 1,     // 304_K_C:305_E_N:305_E_CA:305_E_C [ 304_K_C:305_E_N:305_E_CA -- 305_E_N:305_E_CA:305_E_C ] 
        [ [ -0.6898683669239155, -0.7209624033060767, 0.06553510003697005, 8.250114793224563 ], [ 0.7206393380769293, -0.6925361944465854, -0.032750019808945534, -8.363766010692183 ], [ 0.06899696177205851, 0.024633968429124298, 0.9973126825954125, 21.63972862368898 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 304, ' ') K sidechain

     ],
     [ // 45 : (' ', 305, ' ') E backbone
      [ 144.04292, 180, 181, 0, 0, 1,     // 305_E_N:305_E_CA:305_E_C:305_E_O [ 305_E_N:305_E_CA:305_E_C -- 305_E_CA:305_E_C:305_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -31.22943, 180, 182, 0, 0, 1,     // 305_E_N:305_E_CA:305_E_C:306_L_N [ 305_E_N:305_E_CA:305_E_C -- 305_E_CA:305_E_C:306_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.95105, 182, 183, 0, 0, 1,     // 305_E_CA:305_E_C:306_L_N:306_L_CA [ 305_E_CA:305_E_C:306_L_N -- 305_E_C:306_L_N:306_L_CA ] 
        [ [ 0.4234626422771547, 0.5184663385512889, 0.7428809099612608, 0.0 ], [ -0.2567555013789736, 0.8550980387003705, -0.45042641654593435, 0.0 ], [ -0.8687669440691689, 0.0, 0.4952211595769284, 15.36851666879452 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.66512, 183, 184, 0, 0, 1,     // 305_E_C:306_L_N:306_L_CA:306_L_C [ 305_E_C:306_L_N:306_L_CA -- 306_L_N:306_L_CA:306_L_C ] 
        [ [ -0.8143499318768564, -0.5665137929832987, 0.1260805726979197, 9.949559440715046 ], [ 0.5749331931116411, -0.8171265788061846, 0.041904387204117884, -6.032655227775 ], [ 0.07928437368492992, 0.1066127411155742, 0.9911345577272603, 22.00111726121979 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 305, ' ') E sidechain

     ],
     [ // 46 : (' ', 306, ' ') L backbone
      [ 129.84713, 184, 185, 0, 0, 1,     // 306_L_N:306_L_CA:306_L_C:306_L_O [ 306_L_N:306_L_CA:306_L_C -- 306_L_CA:306_L_C:306_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.27333, 184, 186, 0, 0, 1,     // 306_L_N:306_L_CA:306_L_C:307_V_N [ 306_L_N:306_L_CA:306_L_C -- 306_L_CA:306_L_C:307_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.40422, 186, 187, 0, 0, 1,     // 306_L_CA:306_L_C:307_V_N:307_V_CA [ 306_L_CA:306_L_C:307_V_N -- 306_L_C:307_V_N:307_V_CA ] 
        [ [ 0.33308038488967284, 0.7226455285878349, 0.60567391987253, 0.0 ], [ -0.34822410522726327, 0.6912188076231784, -0.6332112858490941, 0.0 ], [ -0.8762405090729481, 0.0, 0.48187401907509053, 15.356059219125656 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.82803, 187, 188, 0, 0, 1,     // 306_L_C:307_V_N:307_V_CA:307_V_C [ 306_L_C:307_V_N:307_V_CA -- 307_V_N:307_V_CA:307_V_C ] 
        [ [ -0.6575306155289978, -0.7470102986962349, 0.09812799439415436, 8.103546649267459 ], [ 0.7502844458732659, -0.6610947592932999, -0.005193218229600026, -8.471979765614952 ], [ 0.06875129033482599, 0.07020920791958773, 0.9951602520200455, 21.803238905295427 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 306, ' ') L sidechain

     ],
     [ // 47 : (' ', 307, ' ') V backbone
      [ 135.71980, 188, 189, 0, 0, 1,     // 307_V_N:307_V_CA:307_V_C:307_V_O [ 307_V_N:307_V_CA:307_V_C -- 307_V_CA:307_V_C:307_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.45825, 188, 190, 0, 0, 1,     // 307_V_N:307_V_CA:307_V_C:308_H_N [ 307_V_N:307_V_CA:307_V_C -- 307_V_CA:307_V_C:308_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.03555, 190, 191, 0, 0, 1,     // 307_V_CA:307_V_C:308_H_N:308_H_CA [ 307_V_CA:307_V_C:308_H_N -- 307_V_C:308_H_N:308_H_CA ] 
        [ [ 0.3396507948843402, 0.662074073273986, 0.6680533354701768, 0.0 ], [ -0.30005666114507906, 0.749438404072271, -0.5901763114562246, 0.0 ], [ -0.891405260045566, 0.0, 0.45320708551731287, 15.361296004118552 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.35389, 191, 192, 0, 0, 1,     // 307_V_C:308_H_N:308_H_CA:308_H_C [ 307_V_C:308_H_N:308_H_CA -- 308_H_N:308_H_CA:308_H_C ] 
        [ [ -0.7130610008709196, -0.6889829409854353, 0.12979412955917027, 8.911312563695333 ], [ 0.6930549325907015, -0.7206607439239994, -0.01797088140871636, -7.872493556781303 ], [ 0.1059191646901575, 0.07714012702846765, 0.9913780970716325, 21.40672648203654 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 307, ' ') V sidechain

     ],
     [ // 48 : (' ', 308, ' ') H backbone
      [ 139.46450, 192, 193, 0, 0, 1,     // 308_H_N:308_H_CA:308_H_C:308_H_O [ 308_H_N:308_H_CA:308_H_C -- 308_H_CA:308_H_C:308_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -36.26331, 192, 194, 0, 0, 1,     // 308_H_N:308_H_CA:308_H_C:309_M_N [ 308_H_N:308_H_CA:308_H_C -- 308_H_CA:308_H_C:309_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.83308, 194, 195, 0, 0, 1,     // 308_H_CA:308_H_C:309_M_N:309_M_CA [ 308_H_CA:308_H_C:309_M_N -- 308_H_C:309_M_N:309_M_CA ] 
        [ [ 0.3974179562915476, 0.5914969493908455, 0.7015627747239537, 0.0 ], [ -0.29154086467010376, 0.8063072360219294, -0.5146577167070214, 0.0 ], [ -0.8700936111962079, 0.0, 0.49288650595805744, 15.393051747627473 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.55861, 195, 196, 0, 0, 1,     // 308_H_C:309_M_N:309_M_CA:309_M_C [ 308_H_C:309_M_N:309_M_CA -- 309_M_N:309_M_CA:309_M_C ] 
        [ [ -0.7710604329267259, -0.6248839531983933, 0.12241672193780667, 9.375932962770023 ], [ 0.6297756225753367, -0.7767753729553355, 0.0016386518685988832, -6.878067686125426 ], [ 0.09406632758164336, 0.07835856689123212, 0.9924774360206157, 21.980161275973018 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 308, ' ') H sidechain

     ],
     [ // 49 : (' ', 309, ' ') M backbone
      [ 132.80834, 196, 197, 0, 0, 1,     // 309_M_N:309_M_CA:309_M_C:309_M_O [ 309_M_N:309_M_CA:309_M_C -- 309_M_CA:309_M_C:309_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.03659, 196, 198, 0, 0, 1,     // 309_M_N:309_M_CA:309_M_C:310_I_N [ 309_M_N:309_M_CA:309_M_C -- 309_M_CA:309_M_C:310_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.30766, 198, 199, 0, 0, 1,     // 309_M_CA:309_M_C:310_I_N:310_I_CA [ 309_M_CA:309_M_C:310_I_N -- 309_M_C:310_I_N:310_I_CA ] 
        [ [ 0.3389515497741632, 0.7075582502236301, 0.6200590048101615, 0.0 ], [ -0.3393847880798492, 0.7066550237141705, -0.6208515467320034, 0.0 ], [ -0.8774564448026403, 0.0, 0.4796563222499117, 15.328744250112315 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.96799, 199, 200, 0, 0, 1,     // 309_M_C:310_I_N:310_I_CA:310_I_C [ 309_M_C:310_I_N:310_I_CA -- 310_I_N:310_I_CA:310_I_C ] 
        [ [ -0.6790413983485925, -0.7279176516586251, 0.09507192925664616, 8.285070792320528 ], [ 0.7298611415200722, -0.6833321441707636, -0.018970894602093728, -8.295660536002021 ], [ 0.07877495431801453, 0.056507284015241833, 0.9952896228864317, 21.737790090905015 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 309, ' ') M sidechain

     ],
     [ // 50 : (' ', 310, ' ') I backbone
      [ 139.89195, 200, 201, 0, 0, 1,     // 310_I_N:310_I_CA:310_I_C:310_I_O [ 310_I_N:310_I_CA:310_I_C -- 310_I_CA:310_I_C:310_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.56593, 200, 202, 0, 0, 1,     // 310_I_N:310_I_CA:310_I_C:311_S_N [ 310_I_N:310_I_CA:310_I_C -- 310_I_CA:310_I_C:311_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.75133, 202, 203, 0, 0, 1,     // 310_I_CA:310_I_C:311_S_N:311_S_CA [ 310_I_CA:310_I_C:311_S_N -- 310_I_C:311_S_N:311_S_CA ] 
        [ [ 0.35243398582565244, 0.6234148201524048, 0.6979571961441394, 0.0 ], [ -0.28100143538641875, 0.781891272501711, -0.5564927953679433, 0.0 ], [ -0.8926524961852828, 0.0, 0.4507455169540612, 15.36316115718608 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.63942, 203, 204, 0, 0, 1,     // 310_I_C:311_S_N:311_S_CA:311_S_C [ 310_I_C:311_S_N:311_S_CA -- 311_S_N:311_S_CA:311_S_C ] 
        [ [ -0.7563642076100313, -0.6423852779733211, 0.12350846161129339, 9.318116089355602 ], [ 0.6432433018025747, -0.7647105235280277, -0.038155863135456165, -7.4294877949186215 ], [ 0.11895898508549549, 0.0505862614612788, 0.9916096964122527, 21.380864072118356 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 310, ' ') I sidechain

     ],
     [ // 51 : (' ', 311, ' ') S backbone
      [ 138.26210, 204, 205, 0, 0, 1,     // 311_S_N:311_S_CA:311_S_C:311_S_O [ 311_S_N:311_S_CA:311_S_C -- 311_S_CA:311_S_C:311_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.47924, 204, 206, 0, 0, 1,     // 311_S_N:311_S_CA:311_S_C:312_W_N [ 311_S_N:311_S_CA:311_S_C -- 311_S_CA:311_S_C:312_W_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.32613, 206, 207, 0, 0, 1,     // 311_S_CA:311_S_C:312_W_N:312_W_CA [ 311_S_CA:311_S_C:312_W_N -- 311_S_C:312_W_N:312_W_CA ] 
        [ [ 0.3660738455858347, 0.6357986339635621, 0.6795219176951299, 0.0 ], [ -0.3015453155621712, 0.7718549715134757, -0.5597413023999946, 0.0 ], [ -0.8803751258642586, 0.0, 0.47427801736902264, 15.356067215639335 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.64832, 207, 208, 0, 0, 1,     // 311_S_C:312_W_N:312_W_CA:312_W_C [ 311_S_C:312_W_N:312_W_CA -- 312_W_N:312_W_CA:312_W_C ] 
        [ [ -0.7273696371959004, -0.6688759374120662, 0.15342226447499033, 9.078928920930927 ], [ 0.6741388233490154, -0.7382609119568427, -0.022531594046044854, -7.478568926571177 ], [ 0.12833650197472665, 0.08703910746199447, 0.9879038090994033, 21.692781744375516 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 311, ' ') S sidechain

     ],
     [ // 52 : (' ', 312, ' ') W backbone
      [ 126.60081, 208, 209, 0, 0, 1,     // 312_W_N:312_W_CA:312_W_C:312_W_O [ 312_W_N:312_W_CA:312_W_C -- 312_W_CA:312_W_C:312_W_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.23451, 208, 210, 0, 0, 1,     // 312_W_N:312_W_CA:312_W_C:313_A_N [ 312_W_N:312_W_CA:312_W_C -- 312_W_CA:312_W_C:313_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.85708, 210, 211, 0, 0, 1,     // 312_W_CA:312_W_C:313_A_N:313_A_CA [ 312_W_CA:312_W_C:313_A_N -- 312_W_C:313_A_N:313_A_CA ] 
        [ [ 0.31726402888553057, 0.7458773497262411, 0.5856709956457464, 0.0 ], [ -0.355270952739024, 0.6660833124807727, -0.655831968551814, 0.0 ], [ -0.8792758873728015, 0.0, 0.4763128319547696, 15.325229278231733 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -56.51177, 211, 212, 0, 0, 1,     // 312_W_C:313_A_N:313_A_CA:313_A_C [ 312_W_C:313_A_N:313_A_CA -- 313_A_N:313_A_CA:313_A_C ] 
        [ [ -0.6549535180249199, -0.7520572293954229, 0.07379575150956551, 7.829939251470734 ], [ 0.7500431419256666, -0.6588644238626867, -0.057731760915152495, -8.767933722364562 ], [ 0.09203898346383241, 0.017538377409847573, 0.9956009395539832, 21.693139874162437 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 312, ' ') W sidechain

     ],
     [ // 53 : (' ', 313, ' ') A backbone
      [ 133.18132, 212, 213, 0, 0, 1,     // 313_A_N:313_A_CA:313_A_C:313_A_O [ 313_A_N:313_A_CA:313_A_C -- 313_A_CA:313_A_C:313_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.27561, 212, 214, 0, 0, 1,     // 313_A_N:313_A_CA:313_A_C:314_K_N [ 313_A_N:313_A_CA:313_A_C -- 313_A_CA:313_A_C:314_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.76204, 214, 215, 0, 0, 1,     // 313_A_CA:313_A_C:314_K_N:314_K_CA [ 313_A_CA:313_A_C:314_K_N -- 313_A_C:314_K_N:314_K_CA ] 
        [ [ 0.3511467180473296, 0.6596818028008684, 0.6644666293035228, 0.0 ], [ -0.308225180716118, 0.7515450213083685, -0.5832472193839592, 0.0 ], [ -0.8841341642404198, 0.0, 0.4672331105806762, 15.405518264725915 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.73483, 215, 216, 0, 0, 1,     // 313_A_C:314_K_N:314_K_CA:314_K_C [ 313_A_C:314_K_N:314_K_CA -- 314_K_N:314_K_CA:314_K_C ] 
        [ [ -0.7192703044118639, -0.6838274380476996, 0.1225979778152319, 8.9153318866706 ], [ 0.686443160903354, -0.726712572407357, -0.026165319682624335, -7.825588680406897 ], [ 0.10698605535431532, 0.06533660595869212, 0.9921114412612738, 21.674513674825597 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 313, ' ') A sidechain

     ],
     [ // 54 : (' ', 314, ' ') K backbone
      [ 163.53650, 216, 217, 0, 0, 1,     // 314_K_N:314_K_CA:314_K_C:314_K_O [ 314_K_N:314_K_CA:314_K_C -- 314_K_CA:314_K_C:314_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -17.35723, 216, 218, 0, 0, 1,     // 314_K_N:314_K_CA:314_K_C:315_K_N [ 314_K_N:314_K_CA:314_K_C -- 314_K_CA:314_K_C:315_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.76207, 218, 219, 0, 0, 1,     // 314_K_CA:314_K_C:315_K_N:315_K_CA [ 314_K_CA:314_K_C:315_K_N -- 314_K_C:315_K_N:315_K_CA ] 
        [ [ 0.4663767810955983, 0.2983283734057534, 0.8327621987554344, 0.0 ], [ -0.14577137457065595, 0.9544632950622446, -0.2602893096819463, 0.0 ], [ -0.8724926385997133, 0.0, 0.4886272562898125, 15.36011718250815 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -88.94559, 219, 220, 0, 0, 1,     // 314_K_C:315_K_N:315_K_CA:315_K_C [ 314_K_C:315_K_N:315_K_CA -- 315_K_N:315_K_CA:315_K_C ] 
        [ [ -0.9317858109746114, -0.34723733246840166, 0.10583684333262347, 11.133709618253713 ], [ 0.3546605744985354, -0.932973026932567, 0.06145899374910357, -3.47996774477222 ], [ 0.0774020630393947, 0.09480277399165177, 0.9924824203383823, 21.89287526305654 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 314, ' ') K sidechain

     ],
     [ // 55 : (' ', 315, ' ') K backbone
      [ 166.58408, 220, 221, 0, 0, 1,     // 315_K_N:315_K_CA:315_K_C:315_K_O [ 315_K_N:315_K_CA:315_K_C -- 315_K_CA:315_K_C:315_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -14.77122, 220, 222, 0, 0, 1,     // 315_K_N:315_K_CA:315_K_C:316_I_N [ 315_K_N:315_K_CA:315_K_C -- 315_K_CA:315_K_C:316_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.55210, 222, 223, 0, 0, 1,     // 315_K_CA:315_K_C:316_I_N:316_I_CA [ 315_K_CA:315_K_C:316_I_N -- 315_K_C:316_I_N:316_I_CA ] 
        [ [ 0.47090991111264763, 0.2549600800085124, 0.8445349094134187, 0.0 ], [ -0.12416674332728572, 0.9669515797608756, -0.22268197558276462, 0.0 ], [ -0.8733993791315484, 0.0, 0.4870046452885493, 15.336294795993837 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.37085, 223, 224, 0, 0, 1,     // 315_K_C:316_I_N:316_I_CA:316_I_C [ 315_K_C:316_I_N:316_I_CA -- 316_I_N:316_I_CA:316_I_C ] 
        [ [ -0.9621096861201188, -0.26677755649438156, 0.056344362847878554, 11.28885286634279 ], [ 0.26759015607531705, -0.963505421846666, 0.007267079445383051, -2.9765780316723354 ], [ 0.052349405397137005, 0.022068924372635572, 0.9983849469676519, 21.846059785704448 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 315, ' ') K sidechain

     ],
     [ // 56 : (' ', 316, ' ') I backbone
      [ -55.59846, 224, 225, 0, 0, 1,     // 316_I_N:316_I_CA:316_I_C:316_I_O [ 316_I_N:316_I_CA:316_I_C -- 316_I_CA:316_I_C:316_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 124.85201, 224, 226, 0, 0, 1,     // 316_I_N:316_I_CA:316_I_C:317_P_N [ 316_I_N:316_I_CA:316_I_C -- 316_I_CA:316_I_C:317_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.96868, 226, 227, 0, 0, 1,     // 316_I_CA:316_I_C:317_P_N:317_P_CA [ 316_I_CA:316_I_C:317_P_N -- 316_I_C:317_P_N:317_P_CA ] 
        [ [ -0.2784034496059178, -0.8206307672860297, -0.49905577146353836, 0.0 ], [ 0.3997951246385911, -0.5714587857260943, 0.7166580178381257, 0.0 ], [ -0.8733014242303391, 0.0, 0.48718027714313455, 15.372471420784825 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -52.13432, 227, 228, 0, 0, 1,     // 316_I_C:317_P_N:317_P_CA:317_P_C [ 316_I_C:317_P_N:317_P_CA -- 317_P_N:317_P_CA:317_P_C ] 
        [ [ 0.5625793551202544, 0.8255088124142917, -0.04516270384757993, -6.735504979357312 ], [ -0.8247710114761377, 0.5641703003335397, 0.03827075711548199, 9.672373156789968 ], [ 0.057072303450161226, 0.015718551075374017, 0.9982463019370438, 21.947698820558525 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 316, ' ') I sidechain

     ],
     [ // 57 : (' ', 317, ' ') P backbone
      [ -43.05762, 228, 229, 0, 0, 1,     // 317_P_N:317_P_CA:317_P_C:317_P_O [ 317_P_N:317_P_CA:317_P_C -- 317_P_CA:317_P_C:317_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 136.83120, 228, 230, 0, 0, 1,     // 317_P_N:317_P_CA:317_P_C:318_G_N [ 317_P_N:317_P_CA:317_P_C -- 317_P_CA:317_P_C:318_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.30603, 230, 231, 0, 0, 1,     // 317_P_CA:317_P_C:318_G_N:318_G_CA [ 317_P_CA:317_P_C:318_G_N -- 317_P_C:318_G_N:318_G_CA ] 
        [ [ -0.33217623976951716, -0.684150104326767, -0.6493054600742688, 0.0 ], [ 0.3115940778929415, -0.7293412334083914, 0.6090734732956063, 0.0 ], [ -0.8902629254072258, 0.0, 0.45544694932052, 15.311881611897086 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  81.50739, 231, 232, 0, 0, 1,     // 317_P_C:318_G_N:318_G_CA:318_G_C [ 317_P_C:318_G_N:318_G_CA -- 318_G_N:318_G_CA:318_G_C ] 
        [ [ 0.707080366179526, 0.6990067034719031, -0.10689707322827086, -8.661691241808978 ], [ -0.6987154040231406, 0.7138898737789151, 0.04645462621423256, 8.124999239432382 ], [ 0.10878483324529588, 0.04184347759528995, 0.9931842645947053, 21.387513336273678 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 317, ' ') P sidechain

     ],
     [ // 58 : (' ', 318, ' ') G backbone
      [ 164.81508, 232, 233, 0, 0, 1,     // 318_G_N:318_G_CA:318_G_C:318_G_O [ 318_G_N:318_G_CA:318_G_C -- 318_G_CA:318_G_C:318_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -14.64728, 232, 234, 0, 0, 1,     // 318_G_N:318_G_CA:318_G_C:319_F_N [ 318_G_N:318_G_CA:318_G_C -- 318_G_CA:318_G_C:319_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.56745, 234, 235, 0, 0, 1,     // 318_G_CA:318_G_C:319_F_N:319_F_CA [ 318_G_CA:318_G_C:319_F_N -- 318_G_C:319_F_N:319_F_CA ] 
        [ [ 0.478291835952876, 0.25286781705685013, 0.8410081966055566, 0.0 ], [ -0.1250072434148675, 0.9675008357086331, -0.21980746574425694, 0.0 ], [ -0.869258367089235, 0.0, 0.4943580597556358, 15.315075373401092 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -68.59221, 235, 236, 0, 0, 1,     // 318_G_C:319_F_N:319_F_CA:319_F_C [ 318_G_C:319_F_N:319_F_CA -- 319_F_N:319_F_CA:319_F_C ] 
        [ [ -0.9674912216461303, -0.2492498159600014, 0.042839996283655674, 11.264672565373159 ], [ 0.24870035641969365, -0.9684169891788271, -0.01779516188554688, -2.9441557633175295 ], [ 0.04592242104240096, -0.006562340567274892, 0.9989234539902863, 21.93662995001153 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 318, ' ') G sidechain

     ],
     [ // 59 : (' ', 319, ' ') F backbone
      [ 136.33937, 236, 237, 0, 0, 1,     // 319_F_N:319_F_CA:319_F_C:319_F_O [ 319_F_N:319_F_CA:319_F_C -- 319_F_CA:319_F_C:319_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.42287, 236, 238, 0, 0, 1,     // 319_F_N:319_F_CA:319_F_C:320_V_N [ 319_F_N:319_F_CA:319_F_C -- 319_F_CA:319_F_C:320_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.20232, 238, 239, 0, 0, 1,     // 319_F_CA:319_F_C:320_V_N:320_V_CA [ 319_F_CA:319_F_C:320_V_N -- 319_F_C:320_V_N:320_V_CA ] 
        [ [ 0.35239287471376346, 0.6616112379892346, 0.6618835483809233, 0.0 ], [ -0.31092619759344814, 0.7498470309112067, -0.5839985701041808, 0.0 ], [ -0.8826914305129794, 0.0, 0.46995301733146705, 15.359436611632704 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -74.89408, 239, 240, 0, 0, 1,     // 319_F_C:320_V_N:320_V_CA:320_V_C [ 319_F_C:320_V_N:320_V_CA -- 320_V_N:320_V_CA:320_V_C ] 
        [ [ -0.7518459665311469, -0.656641187428924, 0.05958182258691216, 8.86109156918315 ], [ 0.6532148679925159, -0.7541030039996813, -0.06811017245739892, -7.818391647025743 ], [ 0.08965477591498175, -0.012288626065972868, 0.9958970884710152, 21.65102178545158 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 319, ' ') F sidechain

     ],
     [ // 60 : (' ', 320, ' ') V backbone
      [ 163.49866, 240, 241, 0, 0, 1,     // 320_V_N:320_V_CA:320_V_C:320_V_O [ 320_V_N:320_V_CA:320_V_C -- 320_V_CA:320_V_C:320_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -14.81684, 240, 242, 0, 0, 1,     // 320_V_N:320_V_CA:320_V_C:321_E_N [ 320_V_N:320_V_CA:320_V_C -- 320_V_CA:320_V_C:321_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.38212, 242, 243, 0, 0, 1,     // 320_V_CA:320_V_C:321_E_N:321_E_CA [ 320_V_CA:320_V_C:321_E_N -- 320_V_C:321_E_N:321_E_CA ] 
        [ [ 0.44528771025644664, 0.25572986870084047, 0.8580915390264733, 0.0 ], [ -0.11779009096154987, 0.9667482786404384, -0.226987356954124, 0.0 ], [ -0.8876059652603967, 0.0, 0.46060357188601936, 15.34778721997437 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -82.42045, 243, 244, 0, 0, 1,     // 320_V_C:321_E_N:321_E_CA:321_E_C [ 320_V_C:321_E_N:321_E_CA -- 321_E_N:321_E_CA:321_E_C ] 
        [ [ -0.9340331821296768, -0.30534403436405494, 0.18532953180481798, 11.447369183312887 ], [ 0.3186794475820876, -0.9467317199021259, 0.04628671753254961, -3.028124572753906 ], [ 0.16132397332537296, 0.10229404288347113, 0.9815857091569075, 21.49246863766306 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 320, ' ') V sidechain

     ],
     [ // 61 : (' ', 321, ' ') E backbone
      [ 167.82467, 244, 245, 0, 0, 1,     // 321_E_N:321_E_CA:321_E_C:321_E_O [ 321_E_N:321_E_CA:321_E_C -- 321_E_CA:321_E_C:321_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -11.79335, 244, 246, 0, 0, 1,     // 321_E_N:321_E_CA:321_E_C:322_L_N [ 321_E_N:321_E_CA:321_E_C -- 321_E_CA:321_E_C:322_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 167.54610, 246, 247, 0, 0, 1,     // 321_E_CA:321_E_C:322_L_N:322_L_CA [ 321_E_CA:321_E_C:322_L_N -- 321_E_C:322_L_N:322_L_CA ] 
        [ [ 0.4654862891245533, 0.20438241689247139, 0.861133173383911, 0.0 ], [ -0.09718875850062536, 0.9788911214559013, -0.17979576623467833, 0.0 ], [ -0.879702711066733, 0.0, 0.4755240689406164, 15.43554440629864 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.00880, 247, 248, 0, 0, 1,     // 321_E_C:322_L_N:322_L_CA:322_L_C [ 321_E_C:322_L_N:322_L_CA -- 322_L_N:322_L_CA:322_L_C ] 
        [ [ -0.9438522617911709, -0.2999572861472191, 0.13845047634040955, 11.533122595750223 ], [ 0.3197160186393409, -0.9348984846089341, 0.15409896463418274, -2.407997599294595 ], [ 0.08321403329426824, 0.18971149138392707, 0.9783071985321311, 21.804219489121287 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 321, ' ') E sidechain

     ],
     [ // 62 : (' ', 322, ' ') L backbone
      [ -21.37890, 248, 249, 0, 0, 1,     // 322_L_N:322_L_CA:322_L_C:322_L_O [ 322_L_N:322_L_CA:322_L_C -- 322_L_CA:322_L_C:322_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 160.62305, 248, 250, 0, 0, 1,     // 322_L_N:322_L_CA:322_L_C:323_S_N [ 322_L_N:322_L_CA:322_L_C -- 322_L_CA:322_L_C:323_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.46147, 250, 251, 0, 0, 1,     // 322_L_CA:322_L_C:323_S_N:323_S_CA [ 322_L_CA:322_L_C:323_S_N -- 322_L_C:323_S_N:323_S_CA ] 
        [ [ -0.44814947341595723, -0.331781608810043, -0.8301102417946787, 0.0 ], [ 0.15761570189869842, -0.9433562233087879, 0.29195260995504224, 0.0 ], [ -0.8799541692565477, 0.0, 0.4750585858691737, 15.374574542011652 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.80688, 251, 252, 0, 0, 1,     // 322_L_C:323_S_N:323_S_CA:323_S_C [ 322_L_C:323_S_N:323_S_CA -- 323_S_N:323_S_CA:323_S_C ] 
        [ [ 0.9255630643372691, 0.37348587944119743, -0.06197831711687259, -11.10013638373486 ], [ -0.37720325069141053, 0.9237398917343337, -0.06650052696388596, 3.9039559144361 ], [ 0.03241483614700607, 0.0849288542056049, 0.9958596126567669, 21.727002000335116 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 322, ' ') L sidechain

     ],
     [ // 63 : (' ', 323, ' ') S backbone
      [ -29.58849, 252, 253, 0, 0, 1,     // 323_S_N:323_S_CA:323_S_C:323_S_O [ 323_S_N:323_S_CA:323_S_C -- 323_S_CA:323_S_C:323_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 152.13890, 252, 254, 0, 0, 1,     // 323_S_N:323_S_CA:323_S_C:324_L_N [ 323_S_N:323_S_CA:323_S_C -- 323_S_CA:323_S_C:324_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.66687, 254, 255, 0, 0, 1,     // 323_S_CA:323_S_C:324_L_N:324_L_CA [ 323_S_CA:323_S_C:324_L_N -- 323_S_C:324_L_N:324_L_CA ] 
        [ [ -0.41961976924155914, -0.4673296705664159, -0.7781530879396071, 0.0 ], [ 0.22181258988776595, -0.8840831290146223, 0.4113346520280009, 0.0 ], [ -0.8801809042628354, 0.0, 0.4746383631472045, 15.357746459435473 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.70995, 255, 256, 0, 0, 1,     // 323_S_C:324_L_N:324_L_CA:324_L_C [ 323_S_C:324_L_N:324_L_CA -- 324_L_N:324_L_CA:324_L_C ] 
        [ [ 0.8677577302040996, 0.49093634283606885, -0.07731771435948302, -10.419989511404276 ], [ -0.4933679433990281, 0.8696911199881098, -0.015014267882267588, 5.5080456869452705 ], [ 0.05987147983173805, 0.05117482874004448, 0.9968934460644151, 21.713471113563948 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 323, ' ') S sidechain

     ],
     [ // 64 : (' ', 324, ' ') L backbone
      [ 142.78820, 256, 257, 0, 0, 1,     // 324_L_N:324_L_CA:324_L_C:324_L_O [ 324_L_N:324_L_CA:324_L_C -- 324_L_CA:324_L_C:324_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -34.46859, 256, 258, 0, 0, 1,     // 324_L_N:324_L_CA:324_L_C:325_F_N [ 324_L_N:324_L_CA:324_L_C -- 324_L_CA:324_L_C:325_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.66184, 258, 259, 0, 0, 1,     // 324_L_CA:324_L_C:325_F_N:325_F_CA [ 324_L_CA:324_L_C:325_F_N -- 324_L_C:325_F_N:325_F_CA ] 
        [ [ 0.3780432615450245, 0.5659543792807882, 0.7326519862617571, 0.0 ], [ -0.259516923375341, 0.8244365594591847, -0.5029466432027354, 0.0 ], [ -0.888669938099741, 0.0, 0.45854742515665964, 15.37928009307865 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.19634, 259, 260, 0, 0, 1,     // 324_L_C:325_F_N:325_F_CA:325_F_C [ 324_L_C:325_F_N:325_F_CA -- 325_F_N:325_F_CA:325_F_C ] 
        [ [ -0.7836258132065952, -0.604229638253666, 0.14435071573132707, 9.793639741340046 ], [ 0.6120155707001387, -0.7907476115893352, 0.012456162585066943, -6.723080432466156 ], [ 0.1066186010829541, 0.09810585620443826, 0.9894481870626213, 21.50885917088511 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 324, ' ') L sidechain

     ],
     [ // 65 : (' ', 325, ' ') F backbone
      [ 130.83102, 260, 261, 0, 0, 1,     // 325_F_N:325_F_CA:325_F_C:325_F_O [ 325_F_N:325_F_CA:325_F_C -- 325_F_CA:325_F_C:325_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.56883, 260, 262, 0, 0, 1,     // 325_F_N:325_F_CA:325_F_C:326_D_N [ 325_F_N:325_F_CA:325_F_C -- 325_F_CA:325_F_C:326_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.31376, 262, 263, 0, 0, 1,     // 325_F_CA:325_F_C:326_D_N:326_D_CA [ 325_F_CA:325_F_C:326_D_N -- 325_F_C:326_D_N:326_D_CA ] 
        [ [ 0.32716307714777676, 0.7140918835764919, 0.6188999133635226, 0.0 ], [ -0.3337244990048232, 0.7000519850769498, -0.6313122658033834, 0.0 ], [ -0.8840770779265672, 0.0, 0.4673411176911597, 15.325846283441209 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.34033, 263, 264, 0, 0, 1,     // 325_F_C:326_D_N:326_D_CA:326_D_C [ 325_F_C:326_D_N:326_D_CA -- 326_D_N:326_D_CA:326_D_C ] 
        [ [ -0.679342531623022, -0.7286402101022915, 0.0870469353286591, 8.26355332991157 ], [ 0.7290963657363166, -0.6836422300898035, -0.03243132294457662, -8.429283093515853 ], [ 0.083139726954799, 0.04143362716355706, 0.9956761724285426, 21.5657861604937 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 325, ' ') F sidechain

     ],
     [ // 66 : (' ', 326, ' ') D backbone
      [ 130.82696, 264, 265, 0, 0, 1,     // 326_D_N:326_D_CA:326_D_C:326_D_O [ 326_D_N:326_D_CA:326_D_C -- 326_D_CA:326_D_C:326_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.41023, 264, 266, 0, 0, 1,     // 326_D_N:326_D_CA:326_D_C:327_Q_N [ 326_D_N:326_D_CA:326_D_C -- 326_D_CA:326_D_C:327_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.03518, 266, 267, 0, 0, 1,     // 326_D_CA:326_D_C:327_Q_N:327_Q_CA [ 326_D_CA:326_D_C:327_Q_N -- 326_D_C:327_Q_N:327_Q_CA ] 
        [ [ 0.3361721403755545, 0.7121514551246542, 0.6163023584241393, 0.0 ], [ -0.3410208844794765, 0.7020258577601239, -0.625191531760466, 0.0 ], [ -0.877891250887121, 0.0, 0.47886005431216117, 15.32872451256422 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -55.54654, 267, 268, 0, 0, 1,     // 326_D_C:327_Q_N:327_Q_CA:327_Q_C [ 326_D_C:327_Q_N:327_Q_CA -- 327_Q_N:327_Q_CA:327_Q_C ] 
        [ [ -0.6857130649323046, -0.7232587023811564, 0.08181956985386157, 8.249629961237819 ], [ 0.722645936224449, -0.689920926824469, -0.04233161452101133, -8.368617645908547 ], [ 0.08706574205412153, 0.03009923852178997, 0.9957477554084536, 21.738594923241596 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 326, ' ') D sidechain

     ],
     [ // 67 : (' ', 327, ' ') Q backbone
      [ 127.50115, 268, 269, 0, 0, 1,     // 327_Q_N:327_Q_CA:327_Q_C:327_Q_O [ 327_Q_N:327_Q_CA:327_Q_C -- 327_Q_CA:327_Q_C:327_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.77376, 268, 270, 0, 0, 1,     // 327_Q_N:327_Q_CA:327_Q_C:328_V_N [ 327_Q_N:327_Q_CA:327_Q_C -- 327_Q_CA:327_Q_C:328_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.86411, 270, 271, 0, 0, 1,     // 327_Q_CA:327_Q_C:328_V_N:328_V_CA [ 327_Q_CA:327_Q_C:328_V_N -- 327_Q_C:328_V_N:328_V_CA ] 
        [ [ 0.30123914779422606, 0.7635003880194007, 0.5712461232519946, 0.0 ], [ -0.35613747592112066, 0.6458073687209094, -0.6753509759749434, 0.0 ], [ -0.8845456879555412, 0.0, 0.46645356244674385, 15.3659431241576 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.71326, 271, 272, 0, 0, 1,     // 327_Q_C:328_V_N:328_V_CA:328_V_C [ 327_Q_C:328_V_N:328_V_CA -- 328_V_N:328_V_CA:328_V_C ] 
        [ [ -0.6120507410460777, -0.7832381096335861, 0.10923348388911591, 7.654234120331161 ], [ 0.7850128165717818, -0.618440146209501, -0.03587009024813191, -9.04915459920946 ], [ 0.09564919342569099, 0.0637953695339844, 0.9933687042699891, 21.616042035822154 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 327, ' ') Q sidechain

     ],
     [ // 68 : (' ', 328, ' ') V backbone
      [ 127.10746, 272, 273, 0, 0, 1,     // 328_V_N:328_V_CA:328_V_C:328_V_O [ 328_V_N:328_V_CA:328_V_C -- 328_V_CA:328_V_C:328_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.13581, 272, 274, 0, 0, 1,     // 328_V_N:328_V_CA:328_V_C:329_R_N [ 328_V_N:328_V_CA:328_V_C -- 328_V_CA:328_V_C:329_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.60477, 274, 275, 0, 0, 1,     // 328_V_CA:328_V_C:329_R_N:329_R_CA [ 328_V_CA:328_V_C:329_R_N -- 328_V_C:329_R_N:329_R_CA ] 
        [ [ 0.3175499897385254, 0.7447287707454067, 0.5869762022783358, 0.0 ], [ -0.35436052700341936, 0.6673672587294311, -0.6550187469666982, 0.0 ], [ -0.8795400052976108, 0.0, 0.47582494583730967, 15.279568107653398 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.08770, 275, 276, 0, 0, 1,     // 328_V_C:329_R_N:329_R_CA:329_R_C [ 328_V_C:329_R_N:329_R_CA -- 329_R_N:329_R_CA:329_R_C ] 
        [ [ -0.6746079378621695, -0.736775964777683, 0.045445658750484855, 7.856947123566746 ], [ 0.732612067777948, -0.67579770491739, -0.08109882967409413, -8.76772795879367 ], [ 0.09046374035750826, -0.021415876220752254, 0.9956694591711783, 21.648704411658095 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 328, ' ') V sidechain

     ],
     [ // 69 : (' ', 329, ' ') R backbone
      [ 143.91932, 276, 277, 0, 0, 1,     // 329_R_N:329_R_CA:329_R_C:329_R_O [ 329_R_N:329_R_CA:329_R_C -- 329_R_CA:329_R_C:329_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.60293, 276, 278, 0, 0, 1,     // 329_R_N:329_R_CA:329_R_C:330_L_N [ 329_R_N:329_R_CA:329_R_C -- 329_R_CA:329_R_C:330_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.39340, 278, 279, 0, 0, 1,     // 329_R_CA:329_R_C:330_L_N:330_L_CA [ 329_R_CA:329_R_C:330_L_N -- 329_R_C:330_L_N:330_L_CA ] 
        [ [ 0.371199550119804, 0.6101857213024919, 0.6999173376259605, 0.0 ], [ -0.2858924094269447, 0.7922584082977963, -0.539065992911349, 0.0 ], [ -0.8834457675618302, 0.0, 0.46853343079986154, 15.331834766151541 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.06677, 279, 280, 0, 0, 1,     // 329_R_C:330_L_N:330_L_CA:330_L_C [ 329_R_C:330_L_N:330_L_CA -- 330_L_N:330_L_CA:330_L_C ] 
        [ [ -0.7500577959731656, -0.6435319749879378, 0.15257752084760515, 9.365410141390958 ], [ 0.6490725338560089, -0.7605373728754359, -0.016963202924867263, -7.213100527006689 ], [ 0.12695727034564566, 0.08631049546755161, 0.9881459152769536, 21.60115759686944 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 329, ' ') R sidechain

     ],
     [ // 70 : (' ', 330, ' ') L backbone
      [ 121.78602, 280, 281, 0, 0, 1,     // 330_L_N:330_L_CA:330_L_C:330_L_O [ 330_L_N:330_L_CA:330_L_C -- 330_L_CA:330_L_C:330_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -54.26701, 280, 282, 0, 0, 1,     // 330_L_N:330_L_CA:330_L_C:331_L_N [ 330_L_N:330_L_CA:330_L_C -- 330_L_CA:330_L_C:331_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.54143, 282, 283, 0, 0, 1,     // 330_L_CA:330_L_C:331_L_N:331_L_CA [ 330_L_CA:330_L_C:331_L_N -- 330_L_C:331_L_N:331_L_CA ] 
        [ [ 0.27472450783869257, 0.8117473936689334, 0.5153567828839374, 0.0 ], [ -0.38185544301390945, 0.5840087061608702, -0.7163241247989329, 0.0 ], [ -0.8824470893109904, 0.0, 0.4704116649983936, 15.30599232933394 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -55.57628, 283, 284, 0, 0, 1,     // 330_L_C:331_L_N:331_L_CA:331_L_C [ 330_L_C:331_L_N:331_L_CA -- 331_L_N:331_L_CA:331_L_C ] 
        [ [ -0.5691437038148616, -0.8184772387950757, 0.07855223728379554, 6.889709128369396 ], [ 0.8155166214542929, -0.574099707291497, -0.07309012395387555, -9.576404202697551 ], [ 0.10491941926867711, 0.022461871297960814, 0.9942270262863092, 21.594838523021817 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 330, ' ') L sidechain

     ],
     [ // 71 : (' ', 331, ' ') L backbone
      [ 122.25256, 284, 285, 0, 0, 1,     // 331_L_N:331_L_CA:331_L_C:331_L_O [ 331_L_N:331_L_CA:331_L_C -- 331_L_CA:331_L_C:331_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -52.27949, 284, 286, 0, 0, 1,     // 331_L_N:331_L_CA:331_L_C:332_E_N [ 331_L_N:331_L_CA:331_L_C -- 331_L_CA:331_L_C:332_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.59313, 286, 287, 0, 0, 1,     // 331_L_CA:331_L_C:332_E_N:332_E_CA [ 331_L_CA:331_L_C:332_E_N -- 331_L_C:332_E_N:332_E_CA ] 
        [ [ 0.2817990531659821, 0.7910045748228756, 0.5430479318108465, 0.0 ], [ -0.3643357495402927, 0.6118102341488592, -0.7021023422533437, 0.0 ], [ -0.8876084470314989, 0.0, 0.4605987893582992, 15.389850246007564 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -54.34678, 287, 288, 0, 0, 1,     // 331_L_C:332_E_N:332_E_CA:332_E_C [ 331_L_C:332_E_N:332_E_CA -- 332_E_N:332_E_CA:332_E_C ] 
        [ [ -0.5622323950354408, -0.8140384164432208, 0.14572642357282195, 7.2551473785796565 ], [ 0.8159166189469209, -0.5747575623577952, -0.06272013583172517, -9.380122212984075 ], [ 0.13481396403539597, 0.08363731862714475, 0.9873347932864133, 21.543472980646403 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 331, ' ') L sidechain

     ],
     [ // 72 : (' ', 332, ' ') E backbone
      [ 157.55285, 288, 289, 0, 0, 1,     // 332_E_N:332_E_CA:332_E_C:332_E_O [ 332_E_N:332_E_CA:332_E_C -- 332_E_CA:332_E_C:332_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -26.24750, 288, 290, 0, 0, 1,     // 332_E_N:332_E_CA:332_E_C:333_S_N [ 332_E_N:332_E_CA:332_E_C -- 332_E_CA:332_E_C:333_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -174.68029, 290, 291, 0, 0, 1,     // 332_E_CA:332_E_C:333_S_N:333_S_CA [ 332_E_CA:332_E_C:333_S_N -- 332_E_C:333_S_N:333_S_CA ] 
        [ [ 0.44486557667020465, 0.4422494855221182, 0.778787526382715, 0.0 ], [ -0.21935924967309267, 0.8968920740844027, -0.3840131860072464, 0.0 ], [ -0.868317993753869, 0.0, 0.4960079250609368, 15.34874634244328 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -85.07551, 291, 292, 0, 0, 1,     // 332_E_C:333_S_N:333_S_CA:333_S_C [ 332_E_C:333_S_N:333_S_CA -- 333_S_N:333_S_CA:333_S_C ] 
        [ [ -0.9168463798324773, -0.39909981608376216, 0.01058548959604584, 10.426134775535898 ], [ 0.39651841953645195, -0.9133665076000738, -0.09238379599669312, -5.141034103988635 ], [ 0.04653878765497733, -0.08050440731010654, 0.9956672042641825, 21.989126807240122 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 332, ' ') E sidechain

     ],
     [ // 73 : (' ', 333, ' ') S backbone
      [ 156.84797, 292, 293, 0, 0, 1,     // 333_S_N:333_S_CA:333_S_C:333_S_O [ 333_S_N:333_S_CA:333_S_C -- 333_S_CA:333_S_C:333_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -25.85680, 292, 294, 0, 0, 1,     // 333_S_N:333_S_CA:333_S_C:334_C_N [ 333_S_N:333_S_CA:333_S_C -- 333_S_CA:333_S_C:334_C_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -160.37891, 294, 295, 0, 0, 1,     // 333_S_CA:333_S_C:334_C_N:334_C_CA [ 333_S_CA:333_S_C:334_C_N -- 333_S_C:334_C_N:334_C_CA ] 
        [ [ 0.4112626172852393, 0.43612346075513764, 0.8004120105312369, 0.0 ], [ -0.19931536690648802, 0.8998868412078054, -0.38791372434598154, 0.0 ], [ -0.8894585117579273, 0.0, 0.45701592517260614, 15.261488826382008 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -133.84710, 295, 296, 0, 0, 1,     // 333_S_C:334_C_N:334_C_CA:334_C_C [ 333_S_C:334_C_N:334_C_CA -- 334_C_N:334_C_CA:334_C_C ] 
        [ [ -0.9620993014588345, -0.27269819886601276, -0.0007914970975125027, 10.687719120484415 ], [ 0.25947988969508323, -0.914563699074956, -0.3102312479073646, -5.1797235351818465 ], [ 0.0838756280229857, -0.29867864448193454, 0.9506607945814237, 21.363918297163316 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 333, ' ') S sidechain

     ],
     [ // 74 : (' ', 334, ' ') C backbone
      [ -178.15148, 296, 297, 0, 0, 1,     // 334_C_N:334_C_CA:334_C_C:334_C_O [ 334_C_N:334_C_CA:334_C_C -- 334_C_CA:334_C_C:334_C_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   0.54701, 296, 298, 0, 0, 1,     // 334_C_N:334_C_CA:334_C_C:335_W_N [ 334_C_N:334_C_CA:334_C_C -- 334_C_CA:334_C_C:335_W_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -171.65882, 298, 299, 0, 0, 1,     // 334_C_CA:334_C_C:335_W_N:335_W_CA [ 334_C_CA:334_C_C:335_W_N -- 334_C_C:335_W_N:335_W_CA ] 
        [ [ 0.48929414111267977, -0.009546985932407973, 0.8720665677185495, 0.0 ], [ 0.004671497178529559, 0.9999544264913308, 0.008325986698559284, 0.0 ], [ -0.8721063126631501, 0.0, 0.48931644097974475, 15.333886524142235 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -56.48974, 299, 300, 0, 0, 1,     // 334_C_C:335_W_N:335_W_CA:335_W_C [ 334_C_C:335_W_N:335_W_CA -- 335_W_N:335_W_CA:335_W_C ] 
        [ [ -0.9955426972280651, 0.08042661604417446, 0.04925746062833661, 11.64701804182055 ], [ -0.08595101697520835, -0.9886990046086415, -0.12282793235583313, 0.11119898512767534 ], [ 0.03882216733770333, -0.12651417990709427, 0.9912048232356614, 21.869026876053795 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 334, ' ') C sidechain

     ],
     [ // 75 : (' ', 335, ' ') W backbone
      [ 127.87311, 300, 301, 0, 0, 1,     // 335_W_N:335_W_CA:335_W_C:335_W_O [ 335_W_N:335_W_CA:335_W_C -- 335_W_CA:335_W_C:335_W_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.68148, 300, 302, 0, 0, 1,     // 335_W_N:335_W_CA:335_W_C:336_M_N [ 335_W_N:335_W_CA:335_W_C -- 335_W_CA:335_W_C:336_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.01801, 302, 303, 0, 0, 1,     // 335_W_CA:335_W_C:336_M_N:336_M_CA [ 335_W_CA:335_W_C:336_M_N -- 335_W_C:336_M_N:336_M_CA ] 
        [ [ 0.31811449752910775, 0.7275510645631518, 0.6078425905732732, 0.0 ], [ -0.33735638782212185, 0.6860535317692231, -0.6446093539055926, 0.0 ], [ -0.8859987776839274, 0.0, 0.4636875736340005, 15.289328447355583 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.26778, 303, 304, 0, 0, 1,     // 335_W_C:336_M_N:336_M_CA:336_M_C [ 335_W_C:336_M_N:336_M_CA -- 336_M_N:336_M_CA:336_M_C ] 
        [ [ -0.6622552263537193, -0.7431149028955176, 0.09590754017161345, 8.097028467827801 ], [ 0.7435682149925481, -0.6675745735720995, -0.03808540892291052, -8.586795940506718 ], [ 0.09232727018589858, 0.04609153734271664, 0.994661372209359, 21.466078056743 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 335, ' ') W sidechain

     ],
     [ // 76 : (' ', 336, ' ') M backbone
      [ 139.69379, 304, 305, 0, 0, 1,     // 336_M_N:336_M_CA:336_M_C:336_M_O [ 336_M_N:336_M_CA:336_M_C -- 336_M_CA:336_M_C:336_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.14622, 304, 306, 0, 0, 1,     // 336_M_N:336_M_CA:336_M_C:337_E_N [ 336_M_N:336_M_CA:336_M_C -- 336_M_CA:336_M_C:337_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.10726, 306, 307, 0, 0, 1,     // 336_M_CA:336_M_C:337_E_N:337_E_CA [ 336_M_CA:336_M_C:337_E_N -- 336_M_C:337_E_N:337_E_CA ] 
        [ [ 0.3934061543691663, 0.5756651174712457, 0.7168272248117259, 0.0 ], [ -0.27696490138972196, 0.8176855584678115, -0.504659063994064, 0.0 ], [ -0.8766538890021793, 0.0, 0.48112156353395213, 15.398366986591427 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.02763, 307, 308, 0, 0, 1,     // 336_M_C:337_E_N:337_E_CA:337_E_C [ 336_M_C:337_E_N:337_E_CA -- 337_E_N:337_E_CA:337_E_C ] 
        [ [ -0.775746353584658, -0.6187176254456139, 0.12412129093369814, 9.58679849015057 ], [ 0.6266016609400756, -0.7785369219332721, 0.035364101767256434, -6.74927595559698 ], [ 0.07475261471845156, 0.10520818005082753, 0.9916366700778806, 21.83285401389096 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 336, ' ') M sidechain

     ],
     [ // 77 : (' ', 337, ' ') E backbone
      [ 134.45065, 308, 309, 0, 0, 1,     // 337_E_N:337_E_CA:337_E_C:337_E_O [ 337_E_N:337_E_CA:337_E_C -- 337_E_CA:337_E_C:337_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.95397, 308, 310, 0, 0, 1,     // 337_E_N:337_E_CA:337_E_C:338_V_N [ 337_E_N:337_E_CA:337_E_C -- 337_E_CA:337_E_C:338_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.83618, 310, 311, 0, 0, 1,     // 337_E_CA:337_E_C:338_V_N:338_V_CA [ 337_E_CA:337_E_C:338_V_N -- 337_E_C:338_V_N:338_V_CA ] 
        [ [ 0.34232548296502757, 0.6940803005852283, 0.6332975604344925, 0.0 ], [ -0.3300488545510787, 0.7198975873966512, -0.6105859622280894, 0.0 ], [ -0.8797050740573692, 0.0, 0.4755196974655398, 15.319678084462709 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.36320, 311, 312, 0, 0, 1,     // 337_E_C:338_V_N:338_V_CA:338_V_C [ 337_E_C:338_V_N:338_V_CA -- 338_V_N:338_V_CA:338_V_C ] 
        [ [ -0.7091495662612023, -0.7008901197741263, 0.07655019709027133, 8.4697343750949 ], [ 0.6993117418435427, -0.713045428423422, -0.05029219347183715, -8.165989001577051 ], [ 0.08983306958630939, 0.017867764478786356, 0.9957965467911766, 21.67928813079765 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 337, ' ') E sidechain

     ],
     [ // 78 : (' ', 338, ' ') V backbone
      [ 136.34323, 312, 313, 0, 0, 1,     // 338_V_N:338_V_CA:338_V_C:338_V_O [ 338_V_N:338_V_CA:338_V_C -- 338_V_CA:338_V_C:338_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.63854, 312, 314, 0, 0, 1,     // 338_V_N:338_V_CA:338_V_C:339_L_N [ 338_V_N:338_V_CA:338_V_C -- 338_V_CA:338_V_C:339_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.96844, 314, 315, 0, 0, 1,     // 338_V_CA:338_V_C:339_L_N:339_L_CA [ 338_V_CA:338_V_C:339_L_N -- 338_V_C:339_L_N:339_L_CA ] 
        [ [ 0.36392242012893344, 0.6512847736437994, 0.6658743242889343, 0.0 ], [ -0.312344090283331, 0.7588334096623216, -0.5715006786044383, 0.0 ], [ -0.877497374008936, 0.0, 0.4795814410581602, 15.31859477428552 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.42993, 315, 316, 0, 0, 1,     // 338_V_C:339_L_N:339_L_CA:339_L_C [ 338_V_C:339_L_N:339_L_CA -- 339_L_N:339_L_CA:339_L_C ] 
        [ [ -0.7246806415413676, -0.6806927151171885, 0.10721658156079086, 8.919419584472916 ], [ 0.6850123745008793, -0.7285152773985488, 0.004851533528962028, -7.655279921969352 ], [ 0.07480651412717244, 0.07696049755105219, 0.9942238516857437, 21.742611952439663 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 338, ' ') V sidechain

     ],
     [ // 79 : (' ', 339, ' ') L backbone
      [ 127.88520, 316, 317, 0, 0, 1,     // 339_L_N:339_L_CA:339_L_C:339_L_O [ 339_L_N:339_L_CA:339_L_C -- 339_L_CA:339_L_C:339_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.08854, 316, 318, 0, 0, 1,     // 339_L_N:339_L_CA:339_L_C:340_M_N [ 339_L_N:339_L_CA:339_L_C -- 339_L_CA:339_L_C:340_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.89321, 318, 319, 0, 0, 1,     // 339_L_CA:339_L_C:340_M_N:340_M_CA [ 339_L_CA:339_L_C:340_M_N -- 339_L_C:340_M_N:340_M_CA ] 
        [ [ 0.31765820273929113, 0.7441779483159453, 0.5876159013102976, 0.0 ], [ -0.35389341904188504, 0.6679814228257176, -0.6546451456465361, 0.0 ], [ -0.8796889871945016, 0.0, 0.47554945674315713, 15.352338463215661 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.52836, 319, 320, 0, 0, 1,     // 339_L_C:340_M_N:340_M_CA:340_M_C [ 339_L_C:340_M_N:340_M_CA -- 340_M_N:340_M_CA:340_M_C ] 
        [ [ -0.6439407994267394, -0.7603003850366827, 0.0853438418793284, 7.866149698687267 ], [ 0.7613200363090578, -0.6478196175822132, -0.026861596937078985, -8.763440035730598 ], [ 0.07571029750322474, 0.0476766985927833, 0.9959893991721328, 21.718305087893054 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 339, ' ') L sidechain

     ],
     [ // 80 : (' ', 340, ' ') M backbone
      [ 138.22663, 320, 321, 0, 0, 1,     // 340_M_N:340_M_CA:340_M_C:340_M_O [ 340_M_N:340_M_CA:340_M_C -- 340_M_CA:340_M_C:340_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.07828, 320, 322, 0, 0, 1,     // 340_M_N:340_M_CA:340_M_C:341_M_N [ 340_M_N:340_M_CA:340_M_C -- 340_M_CA:340_M_C:341_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.79660, 322, 323, 0, 0, 1,     // 340_M_CA:340_M_C:341_M_N:341_M_CA [ 340_M_CA:340_M_C:341_M_N -- 340_M_C:341_M_N:341_M_CA ] 
        [ [ 0.3695857975583927, 0.6303816052267789, 0.6826605086240424, 0.0 ], [ -0.30012169041182196, 0.776285406143713, -0.554353622837121, 0.0 ], [ -0.8793937168228332, 0.0, 0.4760952539277437, 15.333823148870875 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.16307, 323, 324, 0, 0, 1,     // 340_M_C:341_M_N:341_M_CA:341_M_C [ 340_M_C:341_M_N:341_M_CA -- 341_M_N:341_M_CA:341_M_C ] 
        [ [ -0.7342802101178104, -0.6666272664333319, 0.12822113974289664, 9.116822265159314 ], [ 0.6732618023223508, -0.7393092593519504, 0.01184755544955843, -7.403304259740799 ], [ 0.0868971723533274, 0.09502582114401277, 0.9916748331751201, 21.691999367521213 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 340, ' ') M sidechain

     ],
     [ // 81 : (' ', 341, ' ') M backbone
      [ 132.35150, 324, 325, 0, 0, 1,     // 341_M_N:341_M_CA:341_M_C:341_M_O [ 341_M_N:341_M_CA:341_M_C -- 341_M_CA:341_M_C:341_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.22996, 324, 326, 0, 0, 1,     // 341_M_N:341_M_CA:341_M_C:342_G_N [ 341_M_N:341_M_CA:341_M_C -- 341_M_CA:341_M_C:342_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.90100, 326, 327, 0, 0, 1,     // 341_M_CA:341_M_C:342_G_N:342_G_CA [ 341_M_CA:341_M_C:342_G_N -- 341_M_C:342_G_N:342_G_CA ] 
        [ [ 0.339565515373023, 0.6849281350095707, 0.6446461902794156, 0.0 ], [ -0.3192074966276802, 0.7286106298101278, -0.6059976272431917, 0.0 ], [ -0.8847608913520889, 0.0, 0.4660452393640098, 15.384670215900151 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.12263, 327, 328, 0, 0, 1,     // 341_M_C:342_G_N:342_G_CA:342_G_C [ 341_M_C:342_G_N:342_G_CA -- 342_G_N:342_G_CA:342_G_C ] 
        [ [ -0.6757091865329781, -0.7216424630638562, 0.15049668015631523, 8.633966176709484 ], [ 0.7298787999862805, -0.6835761739781891, -0.0007427650555190605, -8.11633279724509 ], [ 0.10341195562181575, 0.10934244314296052, 0.9886102354123186, 21.626572974592094 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 341, ' ') M sidechain

     ],
     [ // 82 : (' ', 342, ' ') G backbone
      [ 130.94828, 328, 329, 0, 0, 1,     // 342_G_N:342_G_CA:342_G_C:342_G_O [ 342_G_N:342_G_CA:342_G_C -- 342_G_CA:342_G_C:342_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.75019, 328, 330, 0, 0, 1,     // 342_G_N:342_G_CA:342_G_C:343_L_N [ 342_G_N:342_G_CA:342_G_C -- 342_G_CA:342_G_C:343_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.79619, 330, 331, 0, 0, 1,     // 342_G_CA:342_G_C:343_L_N:343_L_CA [ 342_G_CA:342_G_C:343_L_N -- 342_G_C:343_L_N:343_L_CA ] 
        [ [ 0.3271360935067398, 0.7163042586816053, 0.6163523223934064, 0.0 ], [ -0.3358168226196792, 0.697788083156051, -0.6327075569732252, 0.0 ], [ -0.883294423151631, 0.0, 0.468818687798628, 15.32363915845607 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.55362, 331, 332, 0, 0, 1,     // 342_G_C:343_L_N:343_L_CA:343_L_C [ 342_G_C:343_L_N:343_L_CA -- 343_L_N:343_L_CA:343_L_C ] 
        [ [ -0.6700534332067919, -0.7334676828094485, 0.11425215499906138, 8.226862323110996 ], [ 0.7333301925600306, -0.6779293486232578, -0.05136756715198558, -8.445166462256385 ], [ 0.11513133946784278, 0.04936554010019179, 0.9921228845878703, 21.581272429056575 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 342, ' ') G sidechain

     ],
     [ // 83 : (' ', 343, ' ') L backbone
      [ 129.27714, 332, 333, 0, 0, 1,     // 343_L_N:343_L_CA:343_L_C:343_L_O [ 343_L_N:343_L_CA:343_L_C -- 343_L_CA:343_L_C:343_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.19419, 332, 334, 0, 0, 1,     // 343_L_N:343_L_CA:343_L_C:344_M_N [ 343_L_N:343_L_CA:343_L_C -- 343_L_CA:343_L_C:344_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.39390, 334, 335, 0, 0, 1,     // 343_L_CA:343_L_C:344_M_N:344_M_CA [ 343_L_CA:343_L_C:344_M_N -- 343_L_C:344_M_N:344_M_CA ] 
        [ [ 0.3042493577614363, 0.7569288410551138, 0.578352018999429, 0.0 ], [ -0.3524040749333094, 0.6534973064817959, -0.6698900196236963, 0.0 ], [ -0.8850105628025872, 0.0, 0.46557094381828423, 15.34385658945056 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.21900, 335, 336, 0, 0, 1,     // 343_L_C:344_M_N:344_M_CA:344_M_C [ 343_L_C:344_M_N:344_M_CA -- 344_M_N:344_M_CA:344_M_C ] 
        [ [ -0.6552868371940939, -0.7536680686861067, 0.05082915741102703, 7.717099244841407 ], [ 0.7502331423086868, -0.6571885558985595, -0.07248057793367435, -8.938514252113167 ], [ 0.08803063774501861, -0.009361850186779814, 0.9960737736630182, 21.556089037736918 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 343, ' ') L sidechain

     ],
     [ // 84 : (' ', 344, ' ') M backbone
      [ 140.36093, 336, 337, 0, 0, 1,     // 344_M_N:344_M_CA:344_M_C:344_M_O [ 344_M_N:344_M_CA:344_M_C -- 344_M_CA:344_M_C:344_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.70621, 336, 338, 0, 0, 1,     // 344_M_N:344_M_CA:344_M_C:345_W_N [ 344_M_N:344_M_CA:344_M_C -- 344_M_CA:344_M_C:345_W_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.57251, 338, 339, 0, 0, 1,     // 344_M_CA:344_M_C:345_W_N:345_W_CA [ 344_M_CA:344_M_C:345_W_N -- 344_M_C:345_W_N:345_W_CA ] 
        [ [ 0.3828280091530947, 0.5836292267492171, 0.7161142653878562, 0.0 ], [ -0.27515275426439767, 0.8120202741829239, -0.5146980047917469, 0.0 ], [ -0.8818921006725222, 0.0, 0.4714512941666466, 15.320688127379958 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.45700, 339, 340, 0, 0, 1,     // 344_M_C:345_W_N:345_W_CA:345_W_C [ 344_M_C:345_W_N:345_W_CA -- 345_W_N:345_W_CA:345_W_C ] 
        [ [ -0.7725971490641657, -0.6228167037728786, 0.12325988300907142, 9.563173488857208 ], [ 0.6300026588019656, -0.7761139752931457, 0.027271729998748124, -6.873409108148739 ], [ 0.07867842881233425, 0.0987241148664147, 0.9919996239835303, 21.61656957476322 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 344, ' ') M sidechain

     ],
     [ // 85 : (' ', 345, ' ') W backbone
      [ 127.12829, 340, 341, 0, 0, 1,     // 345_W_N:345_W_CA:345_W_C:345_W_O [ 345_W_N:345_W_CA:345_W_C -- 345_W_CA:345_W_C:345_W_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -47.20870, 340, 342, 0, 0, 1,     // 345_W_N:345_W_CA:345_W_C:346_R_N [ 345_W_N:345_W_CA:345_W_C -- 345_W_CA:345_W_C:346_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.28115, 342, 343, 0, 0, 1,     // 345_W_CA:345_W_C:346_R_N:346_R_CA [ 345_W_CA:345_W_C:346_R_N -- 345_W_C:346_R_N:346_R_CA ] 
        [ [ 0.32310335761501874, 0.7338329704340788, 0.5975729175605248, 0.0 ], [ -0.3490261236190615, 0.6793299430349706, -0.6455165323425537, 0.0 ], [ -0.8796504904388748, 0.0, 0.4756206625774863, 15.320295059561989 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.10092, 343, 344, 0, 0, 1,     // 345_W_C:346_R_N:346_R_CA:346_R_C [ 345_W_C:346_R_N:346_R_CA -- 346_R_N:346_R_CA:346_R_C ] 
        [ [ -0.6286599882236746, -0.7665954399481731, 0.13083596851506873, 8.021318994569837 ], [ 0.7734681471177971, -0.6338292709446488, 0.0027350846306277406, -8.664873976093036 ], [ 0.08083096313153149, 0.10291689241539385, 0.9914002565335492, 21.704628992426795 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 345, ' ') W sidechain

     ],
     [ // 86 : (' ', 346, ' ') R backbone
      [ 155.71617, 344, 345, 0, 0, 1,     // 346_R_N:346_R_CA:346_R_C:346_R_O [ 346_R_N:346_R_CA:346_R_C -- 346_R_CA:346_R_C:346_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -25.32447, 344, 346, 0, 0, 1,     // 346_R_N:346_R_CA:346_R_C:347_S_N [ 346_R_N:346_R_CA:346_R_C -- 346_R_CA:346_R_C:347_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -174.98932, 346, 347, 0, 0, 1,     // 346_R_CA:346_R_C:347_S_N:347_S_CA [ 346_R_CA:346_R_C:347_S_N -- 346_R_C:347_S_N:347_S_CA ] 
        [ [ 0.4177632494774042, 0.42774399057716045, 0.801566557380737, 0.0 ], [ -0.19769414075178676, 0.9038999272735484, -0.37931774040675537, 0.0 ], [ -0.8867868368996537, 0.0, 0.462178651499079, 15.325845382001122 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -97.95508, 347, 348, 0, 0, 1,     // 346_R_C:347_S_N:347_S_CA:347_S_C [ 346_R_C:347_S_N:347_S_CA -- 347_S_N:347_S_CA:347_S_C ] 
        [ [ -0.9172698326185315, -0.3896213050321616, 0.08253055696532373, 10.682212314273213 ], [ 0.3791427301878357, -0.9177125025909836, -0.11855189890469109, -5.055042030241153 ], [ 0.12192966953816159, -0.07745321977315951, 0.9895120789728064, 21.48514736008075 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 346, ' ') R sidechain

     ],
     [ // 87 : (' ', 347, ' ') S backbone
      [ -177.49145, 348, 349, 0, 0, 1,     // 347_S_N:347_S_CA:347_S_C:347_S_O [ 347_S_N:347_S_CA:347_S_C -- 347_S_CA:347_S_C:347_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   1.94913, 348, 350, 0, 0, 1,     // 347_S_N:347_S_CA:347_S_C:348_I_N [ 347_S_N:347_S_CA:347_S_C -- 347_S_CA:347_S_C:348_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -170.61813, 350, 351, 0, 0, 1,     // 347_S_CA:347_S_C:348_I_N:348_I_CA [ 347_S_CA:347_S_C:348_I_N -- 347_S_C:348_I_N:348_I_CA ] 
        [ [ 0.4836444925559783, -0.034012196875314225, 0.8746034388704318, 0.0 ], [ 0.0164593347757247, 0.9994214178532072, 0.029764405504322678, 0.0 ], [ -0.8751097617550674, 0.0, 0.48392448262201937, 15.320445669791008 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.80938, 351, 352, 0, 0, 1,     // 347_S_C:348_I_N:348_I_CA:348_I_C [ 347_S_C:348_I_N:348_I_CA -- 348_I_N:348_I_CA:348_I_C ] 
        [ [ -0.9906921926894184, 0.11239792926347886, 0.07678336305160259, 11.72315885132545 ], [ -0.12199467327377671, -0.9833698625227173, -0.1345400058528186, 0.3989612186902607 ], [ 0.06038442710713319, -0.1426548946911181, 0.98792879398447, 21.806955246234445 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 347, ' ') S sidechain

     ],
     [ // 88 : (' ', 348, ' ') I backbone
      [ 127.99579, 352, 353, 0, 0, 1,     // 348_I_N:348_I_CA:348_I_C:348_I_O [ 348_I_N:348_I_CA:348_I_C -- 348_I_CA:348_I_C:348_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.61653, 352, 354, 0, 0, 1,     // 348_I_N:348_I_CA:348_I_C:349_D_N [ 348_I_N:348_I_CA:348_I_C -- 348_I_CA:348_I_C:349_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.55694, 354, 355, 0, 0, 1,     // 348_I_CA:348_I_C:349_D_N:349_D_CA [ 348_I_CA:348_I_C:349_D_N -- 348_I_C:349_D_N:349_D_CA ] 
        [ [ 0.33064412947568306, 0.7023584511437051, 0.6303705781128139, 0.0 ], [ -0.32624761484300285, 0.7118234374527265, -0.6219886556057401, 0.0 ], [ -0.8855715405615288, 0.0, 0.4645030102674046, 15.367429312575094 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -84.92823, 355, 356, 0, 0, 1,     // 348_I_C:349_D_N:349_D_CA:349_D_C [ 348_I_C:349_D_N:349_D_CA -- 349_D_N:349_D_CA:349_D_C ] 
        [ [ -0.7153478804075785, -0.693808926527347, 0.08310585699789838, 8.412564317490443 ], [ 0.6852695101030536, -0.7198137290349019, -0.11078760768240181, -8.300703985420247 ], [ 0.1366861679889586, -0.022301770419388245, 0.9903633285399127, 21.566420168463566 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 348, ' ') I sidechain

     ],
     [ // 89 : (' ', 349, ' ') D backbone
      [ -176.04603, 356, 357, 0, 0, 1,     // 349_D_N:349_D_CA:349_D_C:349_D_O [ 349_D_N:349_D_CA:349_D_C -- 349_D_CA:349_D_C:349_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   6.78934, 356, 358, 0, 0, 1,     // 349_D_N:349_D_CA:349_D_C:350_H_N [ 349_D_N:349_D_CA:349_D_C -- 349_D_CA:349_D_C:350_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -176.99006, 358, 359, 0, 0, 1,     // 349_D_CA:349_D_C:350_H_N:350_H_CA [ 349_D_CA:349_D_C:350_H_N -- 349_D_C:350_H_N:350_H_CA ] 
        [ [ 0.47117975395927675, -0.1182191903497776, 0.8740788651442835, 0.0 ], [ 0.05609585988740542, 0.9929875241074497, 0.10406263244051404, 0.0 ], [ -0.8802516083270555, 0.0, 0.47450722443144316, 15.438952045309996 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -121.55734, 359, 360, 0, 0, 1,     // 349_D_C:350_H_N:350_H_CA:350_H_C [ 349_D_C:350_H_N:350_H_CA -- 350_H_N:350_H_CA:350_H_C ] 
        [ [ -0.9837856437678334, 0.14279741117893413, 0.10851131957959036, 11.70461904586342 ], [ -0.14690170445335746, -0.9886720807172138, -0.030779961647043785, 1.393482348328911 ], [ 0.10288681327073941, -0.046221382182810994, 0.9936185824973799, 21.79298543743792 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 349, ' ') D sidechain

     ],
     [ // 90 : (' ', 350, ' ') H backbone
      [ -104.89876, 360, 361, 0, 0, 1,     // 350_H_N:350_H_CA:350_H_C:350_H_O [ 350_H_N:350_H_CA:350_H_C -- 350_H_CA:350_H_C:350_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  75.41680, 360, 362, 0, 0, 1,     // 350_H_N:350_H_CA:350_H_C:351_P_N [ 350_H_N:350_H_CA:350_H_C -- 350_H_CA:350_H_C:351_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.79758, 362, 363, 0, 0, 1,     // 350_H_CA:350_H_C:351_P_N:351_P_CA [ 350_H_CA:350_H_C:351_P_N -- 350_H_C:351_P_N:351_P_CA ] 
        [ [ 0.12354600718509168, -0.9677830438058432, 0.21939089368184503, 0.0 ], [ 0.4748712343326106, 0.25178558362602366, 0.8432682436102557, 0.0 ], [ -0.8713401717538588, 0.0, 0.4906794321020149, 15.45065457838734 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -53.63930, 363, 364, 0, 0, 1,     // 350_H_C:351_P_N:351_P_CA:351_P_C [ 350_H_C:351_P_N:351_P_CA -- 351_P_N:351_P_CA:351_P_C ] 
        [ [ -0.2714275538140033, 0.9623202512489283, -0.01633453601256556, 2.962375507970389 ], [ -0.9606685267537346, -0.26984892756231454, 0.06555560996800978, 11.386421512749184 ], [ 0.0586776340299731, 0.03348567349879897, 0.9977152123401629, 22.07616536203804 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 350, ' ') H sidechain

     ],
     [ // 91 : (' ', 351, ' ') P backbone
      [ -49.84248, 364, 365, 0, 0, 1,     // 351_P_N:351_P_CA:351_P_C:351_P_O [ 351_P_N:351_P_CA:351_P_C -- 351_P_CA:351_P_C:351_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 131.13603, 364, 366, 0, 0, 1,     // 351_P_N:351_P_CA:351_P_C:352_G_N [ 351_P_N:351_P_CA:351_P_C -- 351_P_CA:351_P_C:352_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.88823, 366, 367, 0, 0, 1,     // 351_P_CA:351_P_C:352_G_N:352_G_CA [ 351_P_CA:351_P_C:352_G_N -- 351_P_C:352_G_N:352_G_CA ] 
        [ [ -0.29448743120434273, -0.7531498501137116, -0.5882537344856905, 0.0 ], [ 0.3371490502771894, -0.6578489973190609, 0.6734725047271756, 0.0 ], [ -0.8942078453915827, 0.0, 0.4476520180230882, 15.373992422569453 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  78.00778, 367, 368, 0, 0, 1,     // 351_P_C:352_G_N:352_G_CA:352_G_C [ 351_P_C:352_G_N:352_G_CA -- 352_G_N:352_G_CA:352_G_C ] 
        [ [ 0.6438969228445296, 0.7587219620018067, -0.09867997328472779, -7.8636175142362434 ], [ -0.7537345982776918, 0.6511835038596232, 0.08856748647370812, 9.002799086621003 ], [ 0.13145686787122338, 0.01735017801732251, 0.9911700475762244, 21.35808405317717 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 351, ' ') P sidechain

     ],
     [ // 92 : (' ', 352, ' ') G backbone
      [ 178.73860, 368, 369, 0, 0, 1,     // 352_G_N:352_G_CA:352_G_C:352_G_O [ 352_G_N:352_G_CA:352_G_C -- 352_G_CA:352_G_C:352_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -0.11589, 368, 370, 0, 0, 1,     // 352_G_N:352_G_CA:352_G_C:353_K_N [ 352_G_N:352_G_CA:352_G_C -- 352_G_CA:352_G_C:353_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.82322, 370, 371, 0, 0, 1,     // 352_G_CA:352_G_C:353_K_N:353_K_CA [ 352_G_CA:352_G_C:353_K_N -- 352_G_C:353_K_N:353_K_CA ] 
        [ [ 0.4873736601002022, 0.0020227195738817556, 0.8731911726798767, 0.0 ], [ -0.0009858222587750065, 0.9999979543006703, -0.0017662244898846056, 0.0 ], [ -0.8731929589701275, 0.0, 0.4873746571222116, 15.360224751916354 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -128.42182, 371, 372, 0, 0, 1,     // 352_G_C:353_K_N:353_K_CA:353_K_C [ 352_G_C:353_K_N:353_K_CA -- 353_K_N:353_K_CA:353_K_C ] 
        [ [ -0.9958058529797554, -0.003526478567367624, 0.09142355889034012, 11.707452004930454 ], [ 0.0037580131395162104, -0.9999901525798902, 0.0023605254692858745, -0.0236809407747387 ], [ 0.09141433426167263, 0.00269419601399389, 0.9958092994139658, 21.894779366303865 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 352, ' ') G sidechain

     ],
     [ // 93 : (' ', 353, ' ') K backbone
      [ -25.22752, 372, 373, 0, 0, 1,     // 353_K_N:353_K_CA:353_K_C:353_K_O [ 353_K_N:353_K_CA:353_K_C -- 353_K_CA:353_K_C:353_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 155.91277, 372, 374, 0, 0, 1,     // 353_K_N:353_K_CA:353_K_C:354_L_N [ 353_K_N:353_K_CA:353_K_C -- 353_K_CA:353_K_C:354_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.66796, 374, 375, 0, 0, 1,     // 353_K_CA:353_K_C:354_L_N:354_L_CA [ 353_K_CA:353_K_C:354_L_N -- 353_K_C:354_L_N:354_L_CA ] 
        [ [ -0.40425568139821655, -0.4081270117395607, -0.818541194043281, 0.0 ], [ 0.1807241940505458, -0.9129251570027723, 0.36593226613164925, 0.0 ], [ -0.8966136903605947, 0.0, 0.4428136066766193, 15.266336594500787 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -119.51870, 375, 376, 0, 0, 1,     // 353_K_C:354_L_N:354_L_CA:354_L_C [ 353_K_C:354_L_N:354_L_CA -- 354_L_N:354_L_CA:354_L_C ] 
        [ [ 0.8800616015654192, 0.4439273996620775, -0.16858244653396, -10.88089945307574 ], [ -0.45160496811362344, 0.8921806191069978, -0.008166741390992164, 4.864351633601588 ], [ 0.146780551249825, 0.08331990589961281, 0.985653723706088, 21.152674833511096 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 353, ' ') K sidechain

     ],
     [ // 94 : (' ', 354, ' ') L backbone
      [ -53.74665, 376, 377, 0, 0, 1,     // 354_L_N:354_L_CA:354_L_C:354_L_O [ 354_L_N:354_L_CA:354_L_C -- 354_L_CA:354_L_C:354_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 125.75499, 376, 378, 0, 0, 1,     // 354_L_N:354_L_CA:354_L_C:355_I_N [ 354_L_N:354_L_CA:354_L_C -- 354_L_CA:354_L_C:355_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -178.51501, 378, 379, 0, 0, 1,     // 354_L_CA:354_L_C:355_I_N:355_I_CA [ 354_L_CA:354_L_C:355_I_N -- 354_L_C:355_I_N:355_I_CA ] 
        [ [ -0.25548967789533933, -0.811523117374137, -0.5255047615924127, 0.0 ], [ 0.354832400497409, -0.584320314525657, 0.7298381584912395, 0.0 ], [ -0.8993436451358191, 0.0, 0.4372425047428691, 15.356295529770177 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -94.68678, 379, 380, 0, 0, 1,     // 354_L_C:355_I_N:355_I_CA:355_I_C [ 354_L_C:355_I_N:355_I_CA -- 355_I_N:355_I_CA:355_I_C ] 
        [ [ 0.5879278285322855, 0.804629516720037, -0.08313969725823614, -7.004324837859883 ], [ -0.7905931510008867, 0.5933195959370671, 0.1514408355350852, 9.727834864231067 ], [ 0.17118217789186094, -0.023306606358608938, 0.9849636866769476, 21.184194266326102 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 354, ' ') L sidechain

     ],
     [ // 95 : (' ', 355, ' ') I backbone
      [ -99.28955, 380, 381, 0, 0, 1,     // 355_I_N:355_I_CA:355_I_C:355_I_O [ 355_I_N:355_I_CA:355_I_C -- 355_I_CA:355_I_C:355_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  79.03572, 380, 382, 0, 0, 1,     // 355_I_N:355_I_CA:355_I_C:356_F_N [ 355_I_N:355_I_CA:355_I_C -- 355_I_CA:355_I_C:356_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.65344, 382, 383, 0, 0, 1,     // 355_I_CA:355_I_C:356_F_N:356_F_CA [ 355_I_CA:355_I_C:356_F_N -- 355_I_C:356_F_N:356_F_CA ] 
        [ [ 0.08240988161048349, -0.9817459464930048, 0.17141618347606574, 0.0 ], [ 0.42537773863042005, 0.19019699404710405, 0.8848044320261489, 0.0 ], [ -0.9012560074088917, 0.0, 0.43328698239086744, 15.333263704242723 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.98814, 383, 384, 0, 0, 1,     // 355_I_C:356_F_N:356_F_CA:356_F_C [ 355_I_C:356_F_N:356_F_CA -- 356_F_N:356_F_CA:356_F_C ] 
        [ [ -0.2105110826894917, 0.9775485469554367, -0.009160917542060696, 2.286379228691885 ], [ -0.9691454302139125, -0.20745403039234897, 0.13311634148165513, 11.801677261829672 ], [ 0.12822721692522324, 0.03690072654142191, 0.9910580796403039, 21.112521744597988 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 355, ' ') I sidechain

     ],
     [ // 96 : (' ', 356, ' ') F backbone
      [ 126.88032, 384, 385, 0, 0, 1,     // 356_F_N:356_F_CA:356_F_C:356_F_O [ 356_F_N:356_F_CA:356_F_C -- 356_F_CA:356_F_C:356_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.51507, 384, 386, 0, 0, 1,     // 356_F_N:356_F_CA:356_F_C:357_A_N [ 356_F_N:356_F_CA:356_F_C -- 356_F_CA:356_F_C:357_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 168.75787, 386, 387, 0, 0, 1,     // 356_F_CA:356_F_C:357_A_N:357_A_CA [ 356_F_CA:356_F_C:357_A_N -- 356_F_C:357_A_N:357_A_CA ] 
        [ [ 0.3175570071102219, 0.7255553655291772, 0.6105136843568869, 0.0 ], [ -0.3348115541131214, 0.6881637970365936, -0.6436861126982113, 0.0 ], [ -0.8871633279546416, 0.0, 0.4614555553164846, 15.27574287609845 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -142.14305, 387, 388, 0, 0, 1,     // 356_F_C:357_A_N:357_A_CA:357_A_C [ 356_F_C:357_A_N:357_A_CA -- 357_A_N:357_A_CA:357_A_C ] 
        [ [ -0.598482074180199, -0.7735429652880806, 0.2084478057890992, 8.106240148066123 ], [ 0.7925371211642345, -0.6096859986571568, 0.012957415565341888, -8.546694927900505 ], [ 0.11706459098153552, 0.17295740485665473, 0.977947655881222, 21.402828386057237 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 356, ' ') F sidechain

     ],
     [ // 97 : (' ', 357, ' ') A backbone
      [ -16.43459, 388, 389, 0, 0, 1,     // 357_A_N:357_A_CA:357_A_C:357_A_O [ 357_A_N:357_A_CA:357_A_C -- 357_A_CA:357_A_C:357_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 165.30196, 388, 390, 0, 0, 1,     // 357_A_N:357_A_CA:357_A_C:358_P_N [ 357_A_N:357_A_CA:357_A_C -- 357_A_CA:357_A_C:358_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.00481, 390, 391, 0, 0, 1,     // 357_A_CA:357_A_C:358_P_N:358_P_CA [ 357_A_CA:357_A_C:358_P_N -- 357_A_C:358_P_N:358_P_CA ] 
        [ [ -0.4885625388701869, -0.2537249289277515, -0.8348235179086203, 0.0 ], [ 0.12815415914194184, -0.9672764136691266, 0.21898139431014316, 0.0 ], [ -0.8630661371571352, 0.0, 0.5050909253715227, 15.425687963245805 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -54.82888, 391, 392, 0, 0, 1,     // 357_A_C:358_P_N:358_P_CA:358_P_C [ 357_A_C:358_P_N:358_P_CA -- 358_P_N:358_P_CA:358_P_C ] 
        [ [ 0.9509723822951774, 0.3033650937204523, -0.06017597547041594, -11.255185891632209 ], [ -0.30737600378814645, 0.9486010142254973, -0.07533995026270834, 2.952332135951225 ], [ 0.03422748029090448, 0.09014286284747791, 0.9953405165422518, 22.235381431722438 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 357, ' ') A sidechain

     ],
     [ // 98 : (' ', 358, ' ') P backbone
      [ 150.37847, 392, 393, 0, 0, 1,     // 358_P_N:358_P_CA:358_P_C:358_P_O [ 358_P_N:358_P_CA:358_P_C -- 358_P_CA:358_P_C:358_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -28.45136, 392, 394, 0, 0, 1,     // 358_P_N:358_P_CA:358_P_C:359_D_N [ 358_P_N:358_P_CA:358_P_C -- 358_P_CA:358_P_C:359_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.89653, 394, 395, 0, 0, 1,     // 358_P_CA:358_P_C:359_D_N:359_D_CA [ 358_P_CA:358_P_C:359_D_N -- 358_P_C:359_D_N:359_D_CA ] 
        [ [ 0.40701775509836857, 0.4764124850587712, 0.7793379825947215, 0.0 ], [ -0.22054539456040254, 0.8792218969521437, -0.4222896930521396, 0.0 ], [ -0.8863951015054634, 0.0, 0.4629295022215793, 15.385367607014262 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -96.60566, 395, 396, 0, 0, 1,     // 358_P_C:359_D_N:359_D_CA:359_D_C [ 358_P_C:359_D_N:359_D_CA -- 359_D_N:359_D_CA:359_D_C ] 
        [ [ -0.8438938118836129, -0.5107299266431485, 0.16431121780216226, 10.424341491941302 ], [ 0.5167307927292503, -0.8561178617049486, -0.007175981824619797, -5.6485017633123515 ], [ 0.14433475710980387, 0.07884889917323387, 0.9863824456007027, 21.577463110669015 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 358, ' ') P sidechain

     ],
     [ // 99 : (' ', 359, ' ') D backbone
      [ 163.22788, 396, 397, 0, 0, 1,     // 359_D_N:359_D_CA:359_D_C:359_D_O [ 359_D_N:359_D_CA:359_D_C -- 359_D_CA:359_D_C:359_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -14.01948, 396, 398, 0, 0, 1,     // 359_D_N:359_D_CA:359_D_C:360_L_N [ 359_D_N:359_D_CA:359_D_C -- 359_D_CA:359_D_C:360_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.99062, 398, 399, 0, 0, 1,     // 359_D_CA:359_D_C:360_L_N:360_L_CA [ 359_D_CA:359_D_C:360_L_N -- 359_D_C:360_L_N:360_L_CA ] 
        [ [ 0.47574880120920865, 0.2422517801055033, 0.8455632165507887, 0.0 ], [ -0.11878932194238678, 0.9702134172622614, -0.21112797530586305, 0.0 ], [ -0.8715229056889262, 0.0, 0.49035479487767947, 15.43282693630014 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -118.22546, 399, 400, 0, 0, 1,     // 359_D_C:360_L_N:360_L_CA:360_L_C [ 359_D_C:360_L_N:360_L_CA -- 360_L_N:360_L_CA:360_L_C ] 
        [ [ -0.9630045123606494, -0.25059501699699077, 0.09911834658278713, 11.322438200551886 ], [ 0.2508704577371536, -0.9679702621944454, -0.009878509096399049, -2.8270901643043445 ], [ 0.0984191170849374, 0.015352816142145636, 0.9950266169448568, 21.99887871418438 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 359, ' ') D sidechain

     ],
     [ // 100 : (' ', 360, ' ') L backbone
      [ -93.76584, 400, 401, 0, 0, 1,     // 360_L_N:360_L_CA:360_L_C:360_L_O [ 360_L_N:360_L_CA:360_L_C -- 360_L_CA:360_L_C:360_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  87.71160, 400, 402, 0, 0, 1,     // 360_L_N:360_L_CA:360_L_C:361_V_N [ 360_L_N:360_L_CA:360_L_C -- 360_L_CA:360_L_C:361_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 167.57920, 402, 403, 0, 0, 1,     // 360_L_CA:360_L_C:361_V_N:361_V_CA [ 360_L_CA:360_L_C:361_V_N -- 360_L_C:361_V_N:361_V_CA ] 
        [ [ 0.016807049778002285, -0.9992024989824095, 0.036220009733123296, 0.0 ], [ 0.4205822905782258, 0.03992951323654964, 0.9063752924835633, 0.0 ], [ -0.9070987046235556, 0.0, 0.42091797309008727, 15.346154331320665 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -97.65203, 403, 404, 0, 0, 1,     // 360_L_C:361_V_N:361_V_CA:361_V_C [ 360_L_C:361_V_N:361_V_CA -- 361_V_N:361_V_CA:361_V_C ] 
        [ [ -0.1610961491776575, 0.9722003841441946, -0.16992481805201354, 0.4819304523832248 ], [ -0.9745132879436001, -0.12945794857788712, 0.1832061439236106, 12.059904399642322 ], [ 0.1561149651426265, 0.19510779743486711, 0.9682773698887274, 20.94673783977914 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 360, ' ') L sidechain

     ],
     [ // 101 : (' ', 361, ' ') V backbone
      [ -52.73159, 404, 405, 0, 0, 1,     // 361_V_N:361_V_CA:361_V_C:361_V_O [ 361_V_N:361_V_CA:361_V_C -- 361_V_CA:361_V_C:361_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 128.16964, 404, 406, 0, 0, 1,     // 361_V_N:361_V_CA:361_V_C:362_L_N [ 361_V_N:361_V_CA:361_V_C -- 361_V_CA:361_V_C:362_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.79702, 406, 407, 0, 0, 1,     // 361_V_CA:361_V_C:362_L_N:362_L_CA [ 361_V_CA:361_V_C:362_L_N -- 361_V_C:362_L_N:362_L_CA ] 
        [ [ -0.2779811133190949, -0.7861844214260116, -0.5519425297482736, 0.0 ], [ 0.3536363527748774, -0.6179919542413537, 0.7021590094053022, 0.0 ], [ -0.893122517146421, 0.0, 0.4498134828637766, 15.178812259721473 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -118.95632, 407, 408, 0, 0, 1,     // 361_V_C:362_L_N:362_L_CA:362_L_C [ 361_V_C:362_L_N:362_L_CA -- 362_L_N:362_L_CA:362_L_C ] 
        [ [ 0.5795603497650816, 0.8044433601744908, -0.1303099430256593, -7.353979222393942 ], [ -0.8042420635225302, 0.5904119428345308, 0.06788549933069014, 9.355435552933953 ], [ 0.13154658582114007, 0.06545699374038759, 0.9891465400683926, 21.172043128235124 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 361, ' ') V sidechain

     ],
     [ // 102 : (' ', 362, ' ') L backbone
      [ -45.33305, 408, 409, 0, 0, 1,     // 362_L_N:362_L_CA:362_L_C:362_L_O [ 362_L_N:362_L_CA:362_L_C -- 362_L_CA:362_L_C:362_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 136.47098, 408, 410, 0, 0, 1,     // 362_L_N:362_L_CA:362_L_C:363_D_N [ 362_L_N:362_L_CA:362_L_C -- 362_L_CA:362_L_C:363_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 161.95108, 410, 411, 0, 0, 1,     // 362_L_CA:362_L_C:363_D_N:363_D_CA [ 362_L_CA:362_L_C:363_D_N -- 362_L_C:363_D_N:363_D_CA ] 
        [ [ -0.3309798714588946, -0.6887218243233553, -0.6450694329990865, 0.0 ], [ 0.31440687490495567, -0.725025688304014, 0.6127691802887856, 0.0 ], [ -0.8897194174016622, 0.0, 0.4565077855836052, 15.27202008119638 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -100.27691, 411, 412, 0, 0, 1,     // 362_L_C:363_D_N:363_D_CA:363_D_C [ 362_L_C:363_D_N:363_D_CA -- 363_D_N:363_D_CA:363_D_C ] 
        [ [ 0.5945210875079624, 0.7573785578366298, -0.27004147206985923, -8.55897243870165 ], [ -0.7997030387774394, 0.5919365212461892, -0.10043009800398048, 8.130403111791757 ], [ 0.08378380677958629, 0.27566079689405715, 0.9575966785537919, 21.329100691413934 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 362, ' ') L sidechain

     ],
     [ // 103 : (' ', 363, ' ') D backbone
      [   0.14976, 412, 413, 0, 0, 1,     // 363_D_N:363_D_CA:363_D_C:363_D_O [ 363_D_N:363_D_CA:363_D_C -- 363_D_CA:363_D_C:363_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.80474, 412, 414, 0, 0, 1,     // 363_D_N:363_D_CA:363_D_C:364_R_N [ 363_D_N:363_D_CA:363_D_C -- 363_D_CA:363_D_C:364_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -176.83085, 414, 415, 0, 0, 1,     // 363_D_CA:363_D_C:364_R_N:364_R_CA [ 363_D_CA:363_D_C:364_R_N -- 363_D_C:364_R_N:364_R_CA ] 
        [ [ -0.43470635751784975, -0.05573885489312446, -0.8988456835234658, 0.0 ], [ 0.024267761691453886, -0.9984453816084299, 0.0501786377583488, 0.0 ], [ -0.9002452213014241, 0.0, 0.4353832122670209, 15.417239666620056 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.04086, 415, 416, 0, 0, 1,     // 363_D_C:364_R_N:364_R_CA:364_R_C [ 363_D_C:364_R_N:364_R_CA -- 364_R_N:364_R_CA:364_R_C ] 
        [ [ 0.9910898807488719, 0.031621380843401777, -0.12938677115819605, -12.025798923260124 ], [ -0.02487051915040908, 0.998260053652309, 0.05346328234665681, 0.6713479510293014 ], [ 0.13085222793050175, -0.04976900195549185, 0.9901518776884576, 21.242300699078786 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 363, ' ') D sidechain

     ],
     [ // 104 : (' ', 364, ' ') R backbone
      [ 155.96791, 416, 417, 0, 0, 1,     // 364_R_N:364_R_CA:364_R_C:364_R_O [ 364_R_N:364_R_CA:364_R_C -- 364_R_CA:364_R_C:364_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -25.17162, 416, 418, 0, 0, 1,     // 364_R_N:364_R_CA:364_R_C:365_D_N [ 364_R_N:364_R_CA:364_R_C -- 364_R_CA:364_R_C:365_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.70794, 418, 419, 0, 0, 1,     // 364_R_CA:364_R_C:365_D_N:365_D_CA [ 364_R_CA:364_R_C:365_D_N -- 364_R_C:365_D_N:365_D_CA ] 
        [ [ 0.4368549090780824, 0.4253310807573207, 0.7926230252498283, 0.0 ], [ -0.20530409302112157, 0.9050378289009855, -0.37250068140469905, 0.0 ], [ -0.8757899393137347, 0.0, 0.48269243022534025, 15.317158773652574 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -80.17992, 419, 420, 0, 0, 1,     // 364_R_C:365_D_N:365_D_CA:365_D_C [ 364_R_C:365_D_N:365_D_CA -- 365_D_N:365_D_CA:365_D_C ] 
        [ [ -0.8971449437541509, -0.4350735185549802, 0.0764328682475726, 10.55987541079545 ], [ 0.435260659621953, -0.9001783419297329, -0.015070199272050725, -4.962713245467737 ], [ 0.07535985723064943, 0.019748067571948486, 0.9969608346095397, 21.747923162949846 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 364, ' ') R sidechain

     ],
     [ // 105 : (' ', 365, ' ') D backbone
      [ 133.09782, 420, 421, 0, 0, 1,     // 365_D_N:365_D_CA:365_D_C:365_D_O [ 365_D_N:365_D_CA:365_D_C -- 365_D_CA:365_D_C:365_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.87246, 420, 422, 0, 0, 1,     // 365_D_N:365_D_CA:365_D_C:366_E_N [ 365_D_N:365_D_CA:365_D_C -- 365_D_CA:365_D_C:366_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 171.40874, 422, 423, 0, 0, 1,     // 365_D_CA:365_D_C:366_E_N:366_E_CA [ 365_D_CA:365_D_C:366_E_N -- 365_D_C:366_E_N:366_E_CA ] 
        [ [ 0.35640816163597094, 0.6543774255062581, 0.666905845908603, 0.0 ], [ -0.3084307012898472, 0.7561680930770629, -0.5771311094665158, 0.0 ], [ -0.8819544913549229, 0.0, 0.47133456819851366, 15.313203468232704 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.47389, 423, 424, 0, 0, 1,     // 365_D_C:366_E_N:366_E_CA:366_E_C [ 365_D_C:366_E_N:366_E_CA -- 366_E_N:366_E_CA:366_E_C ] 
        [ [ -0.6964434734030165, -0.700276626830607, 0.15677733978189787, 8.905677216452002 ], [ 0.7119762478267354, -0.7016085088801363, 0.028901951448903906, -7.706850080882331 ], [ 0.08975695452111143, 0.13175031757737282, 0.9872109921052118, 21.607275639114537 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 365, ' ') D sidechain

     ],
     [ // 106 : (' ', 366, ' ') E backbone
      [ 140.84407, 424, 425, 0, 0, 1,     // 366_E_N:366_E_CA:366_E_C:366_E_O [ 366_E_N:366_E_CA:366_E_C -- 366_E_CA:366_E_C:366_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.49929, 424, 426, 0, 0, 1,     // 366_E_N:366_E_CA:366_E_C:367_G_N [ 366_E_N:366_E_CA:366_E_C -- 366_E_CA:366_E_C:367_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.95263, 426, 427, 0, 0, 1,     // 366_E_CA:366_E_C:367_G_N:367_G_CA [ 366_E_CA:366_E_C:367_G_N -- 366_E_C:367_G_N:367_G_CA ] 
        [ [ 0.3787879033349006, 0.5806928642303265, 0.720635498514422, 0.0 ], [ -0.27017970151207693, 0.8141227164451189, -0.5140108281533234, 0.0 ], [ -0.8851681496630993, 0.0, 0.4652712615475031, 15.400294644666488 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.49498, 427, 428, 0, 0, 1,     // 366_E_C:367_G_N:367_G_CA:367_G_C [ 366_E_C:367_G_N:367_G_CA -- 367_G_N:367_G_CA:367_G_C ] 
        [ [ -0.7801034582764435, -0.6117666004623177, 0.13107334185070596, 9.64192022824887 ], [ 0.6166917893277851, -0.7871955820887693, -0.0037884714091923136, -6.877334535597122 ], [ 0.105498015909437, 0.07787645407119198, 0.9913654354174699, 21.62550613866835 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 366, ' ') E sidechain

     ],
     [ // 107 : (' ', 367, ' ') G backbone
      [ 157.61553, 428, 429, 0, 0, 1,     // 367_G_N:367_G_CA:367_G_C:367_G_O [ 367_G_N:367_G_CA:367_G_C -- 367_G_CA:367_G_C:367_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -19.59498, 428, 430, 0, 0, 1,     // 367_G_N:367_G_CA:367_G_C:368_K_N [ 367_G_N:367_G_CA:367_G_C -- 367_G_CA:367_G_C:368_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.52398, 430, 431, 0, 0, 1,     // 367_G_CA:367_G_C:368_K_N:368_K_CA [ 367_G_CA:367_G_C:368_K_N -- 367_G_C:368_K_N:368_K_CA ] 
        [ [ 0.4705554187001637, 0.3353690227387929, 0.8161526919144118, 0.0 ], [ -0.16751078967952665, 0.9420868423809066, -0.29053832234959015, 0.0 ], [ -0.8663242656608754, 0.0, 0.49948199840148894, 15.376653052413397 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.23094, 431, 432, 0, 0, 1,     // 367_G_C:368_K_N:368_K_CA:368_K_C [ 367_G_C:368_K_N:368_K_CA -- 368_K_N:368_K_CA:368_K_C ] 
        [ [ -0.9481437762046282, -0.3147273715057571, 0.04438537224700485, 10.958389863019217 ], [ 0.3122947855344972, -0.9484439672504894, -0.054092586498658626, -3.9010251856012266 ], [ 0.05912145610850969, -0.03742622892076647, 0.9975489616136042, 22.08314133857524 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 367, ' ') G sidechain

     ],
     [ // 108 : (' ', 368, ' ') K backbone
      [ 165.92037, 432, 433, 0, 0, 1,     // 368_K_N:368_K_CA:368_K_C:368_K_O [ 368_K_N:368_K_CA:368_K_C -- 368_K_CA:368_K_C:368_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -15.44118, 432, 434, 0, 0, 1,     // 368_K_N:368_K_CA:368_K_C:369_C_N [ 368_K_N:368_K_CA:368_K_C -- 368_K_CA:368_K_C:369_C_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.58622, 434, 435, 0, 0, 1,     // 368_K_CA:368_K_C:369_C_N:369_C_CA [ 368_K_CA:368_K_C:369_C_N -- 368_K_C:369_C_N:369_C_CA ] 
        [ [ 0.4610551639882274, 0.26624895633181905, 0.8464866384131565, 0.0 ], [ -0.12735232806810576, 0.9639042967287869, -0.23381593462254974, 0.0 ], [ -0.8781853564569512, 0.0, 0.4783204780318084, 15.306286583668172 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -74.13794, 435, 436, 0, 0, 1,     // 368_K_C:369_C_N:369_C_CA:369_C_C [ 368_K_C:369_C_N:369_C_CA -- 369_C_N:369_C_CA:369_C_C ] 
        [ [ -0.9619720046661979, -0.26291235662065693, 0.07407398311599621, 11.286089210226274 ], [ 0.26162602535192886, -0.9647988757081345, -0.026738591789235693, -3.117435500068156 ], [ 0.07849640180956219, -0.006342094960818223, 0.9968942234432195, 21.683667325631177 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 368, ' ') K sidechain

     ],
     [ // 109 : (' ', 369, ' ') C backbone
      [ 174.50151, 436, 437, 0, 0, 1,     // 369_C_N:369_C_CA:369_C_C:369_C_O [ 369_C_N:369_C_CA:369_C_C -- 369_C_CA:369_C_C:369_C_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -5.71816, 436, 438, 0, 0, 1,     // 369_C_N:369_C_CA:369_C_C:370_V_N [ 369_C_N:369_C_CA:369_C_C -- 369_C_CA:369_C_C:370_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.05788, 438, 439, 0, 0, 1,     // 369_C_CA:369_C_C:370_V_N:370_V_CA [ 369_C_CA:369_C_C:370_V_N -- 369_C_C:370_V_N:370_V_CA ] 
        [ [ 0.49448956095799984, 0.09963505087360144, 0.8634540698502611, 0.0 ], [ -0.049514876196340395, 0.9950240482709023, -0.08646051351828898, 0.0 ], [ -0.8677720617413456, 0.0, 0.49696242198095253, 15.41019254506266 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -118.41976, 439, 440, 0, 0, 1,     // 369_C_C:370_V_N:370_V_CA:370_V_C [ 369_C_C:370_V_N:370_V_CA -- 370_V_N:370_V_CA:370_V_C ] 
        [ [ -0.9894818519663949, -0.11633602637640592, 0.08597437755576558, 11.57585943782862 ], [ 0.11841162121555843, -0.9927744235850792, 0.019432751642957737, -1.1591291145153342 ], [ 0.08309243401330986, 0.029408720513853546, 0.9961078127228407, 22.07269764360486 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 369, ' ') C sidechain

     ],
     [ // 110 : (' ', 370, ' ') V backbone
      [ -53.35889, 440, 441, 0, 0, 1,     // 370_V_N:370_V_CA:370_V_C:370_V_O [ 370_V_N:370_V_CA:370_V_C -- 370_V_CA:370_V_C:370_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 128.63382, 440, 442, 0, 0, 1,     // 370_V_N:370_V_CA:370_V_C:371_E_N [ 370_V_N:370_V_CA:370_V_C -- 370_V_CA:370_V_C:371_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.69665, 442, 443, 0, 0, 1,     // 370_V_CA:370_V_C:371_E_N:371_E_CA [ 370_V_CA:370_V_C:371_E_N -- 370_V_C:371_E_N:371_E_CA ] 
        [ [ -0.28989494056375303, -0.7811520969254085, -0.5529577966758179, 0.0 ], [ 0.3627058336039341, -0.6243407735115792, 0.6918403549957447, 0.0 ], [ -0.8856666425383838, 0.0, 0.4643216539154583, 15.303244410486883 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.00945, 443, 444, 0, 0, 1,     // 370_V_C:371_E_N:371_E_CA:371_E_C [ 370_V_C:371_E_N:371_E_CA -- 371_E_N:371_E_CA:371_E_C ] 
        [ [ 0.6389367750668251, 0.7688700420472717, -0.024467445911498936, -7.376484761232356 ], [ -0.7649852520228078, 0.6384135405541222, 0.08500420827671692, 9.229185059892792 ], [ 0.0809775379645909, -0.0355950794264718, 0.9960801316489635, 21.49731857854783 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 370, ' ') V sidechain

     ],
     [ // 111 : (' ', 371, ' ') E backbone
      [ -56.34496, 444, 445, 0, 0, 1,     // 371_E_N:371_E_CA:371_E_C:371_E_O [ 371_E_N:371_E_CA:371_E_C -- 371_E_CA:371_E_C:371_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 124.96478, 444, 446, 0, 0, 1,     // 371_E_N:371_E_CA:371_E_C:372_G_N [ 371_E_N:371_E_CA:371_E_C -- 371_E_CA:371_E_C:372_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.19251, 446, 447, 0, 0, 1,     // 371_E_CA:371_E_C:372_G_N:372_G_CA [ 371_E_CA:371_E_C:372_G_N -- 371_E_C:372_G_N:372_G_CA ] 
        [ [ -0.25968914098497903, -0.8195044503351798, -0.510856149943718, 0.0 ], [ 0.3713601477573138, -0.5730728190036889, 0.7305334932614871, 0.0 ], [ -0.8914332228003115, 0.0, 0.453152081853157, 15.422375974015793 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  79.13151, 447, 448, 0, 0, 1,     // 371_E_C:372_G_N:372_G_CA:372_G_C [ 371_E_C:372_G_N:372_G_CA -- 372_G_N:372_G_CA:372_G_C ] 
        [ [ 0.5302231460198712, 0.8383854085361319, -0.1263856090640617, -6.827922179031803 ], [ -0.8396207638460788, 0.53993372604884, 0.05923296711750572, 9.764051664475858 ], [ 0.11789990815652819, 0.07470929144838392, 0.9902111559803609, 21.47904603759531 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 371, ' ') E sidechain

     ],
     [ // 112 : (' ', 372, ' ') G backbone
      [ 169.80772, 448, 449, 0, 0, 1,     // 372_G_N:372_G_CA:372_G_C:372_G_O [ 372_G_N:372_G_CA:372_G_C -- 372_G_CA:372_G_C:372_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -6.39587, 448, 450, 0, 0, 1,     // 372_G_N:372_G_CA:372_G_C:373_I_N [ 372_G_N:372_G_CA:372_G_C -- 372_G_CA:372_G_C:373_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.66929, 450, 451, 0, 0, 1,     // 372_G_CA:372_G_C:373_I_N:373_I_CA [ 372_G_CA:372_G_C:373_I_N -- 372_G_C:373_I_N:373_I_CA ] 
        [ [ 0.4999899988688681, 0.11139724239301733, 0.8588368037167134, 0.0 ], [ -0.05604634189275294, 0.993775957843231, -0.0962712479051212, 0.0 ], [ -0.864215718782961, 0.0, 0.5031214479710143, 15.355364788273576 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.66929, 451, 452, 0, 0, 1,     // 372_G_C:373_I_N:373_I_CA:373_I_C [ 372_G_C:373_I_N:373_I_CA -- 373_I_N:373_I_CA:373_I_C ] 
        [ [ -0.9937367414302538, -0.1085094874763368, 0.02670168270359777, 11.520473640375952 ], [ 0.10833419171763836, -0.9940828985549145, -0.00793055518520573, -1.2913866394841644 ], [ 0.027404226616835576, -0.004988178854284035, 0.9996119879409464, 22.104257423258655 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 372, ' ') G sidechain

     ],
     [ // 113 : (' ', 373, ' ') I backbone
      [ 155.69609, 452, 453, 0, 0, 1,     // 373_I_N:373_I_CA:373_I_C:373_I_O [ 373_I_N:373_I_CA:373_I_C -- 373_I_CA:373_I_C:373_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -23.85555, 452, 454, 0, 0, 1,     // 373_I_N:373_I_CA:373_I_C:374_L_N [ 373_I_N:373_I_CA:373_I_C -- 373_I_CA:373_I_C:374_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -173.62701, 454, 455, 0, 0, 1,     // 373_I_CA:373_I_C:374_L_N:374_L_CA [ 373_I_CA:373_I_C:374_L_N -- 373_I_C:374_L_N:374_L_CA ] 
        [ [ 0.4294948890588044, 0.4044321832255869, 0.8074458182712637, 0.0 ], [ -0.18992743839546888, 0.9145679904585472, -0.3570615618813513, 0.0 ], [ -0.8828712864381199, 0.0, 0.46961504616345007, 15.397938581881323 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.15723, 455, 456, 0, 0, 1,     // 373_I_C:374_L_N:374_L_CA:374_L_C [ 373_I_C:374_L_N:374_L_CA -- 374_L_N:374_L_CA:374_L_C ] 
        [ [ -0.9333424831872105, -0.35425878579935893, 0.0580734169986209, 10.814448782729135 ], [ 0.3451202902567228, -0.9299983215832691, -0.1264717640638324, -4.782270074188292 ], [ 0.09881191391247898, -0.09799915579364057, 0.9902688378075353, 21.687683044496683 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 373, ' ') I sidechain

     ],
     [ // 114 : (' ', 374, ' ') L backbone
      [ 143.15369, 456, 457, 0, 0, 1,     // 374_L_N:374_L_CA:374_L_C:374_L_O [ 374_L_N:374_L_CA:374_L_C -- 374_L_CA:374_L_C:374_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.62491, 456, 458, 0, 0, 1,     // 374_L_N:374_L_CA:374_L_C:375_E_N [ 374_L_N:374_L_CA:374_L_C -- 374_L_CA:374_L_C:375_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.22742, 458, 459, 0, 0, 1,     // 374_L_CA:374_L_C:375_E_N:375_E_CA [ 374_L_CA:374_L_C:375_E_N -- 374_L_C:375_E_N:375_E_CA ] 
        [ [ 0.4057536079866497, 0.5537536487757682, 0.7271319729411845, 0.0 ], [ -0.26983642314123985, 0.8326805488706478, -0.4835611732540347, 0.0 ], [ -0.8732424144257755, 0.0, 0.4872860408916328, 15.353067293481368 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.40194, 459, 460, 0, 0, 1,     // 374_L_C:375_E_N:375_E_CA:375_E_C [ 374_L_C:375_E_N:375_E_CA -- 375_E_N:375_E_CA:375_E_C ] 
        [ [ -0.8004367746889645, -0.5917562944777037, 0.09552726140502929, 9.738935806781388 ], [ 0.5970000280585143, -0.801317809394229, 0.03848031764119132, -6.476638904934887 ], [ 0.05377672566883157, 0.08783083908087758, 0.9946827672593349, 21.87959569325526 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 374, ' ') L sidechain

     ],
     [ // 115 : (' ', 375, ' ') E backbone
      [ 127.14783, 460, 461, 0, 0, 1,     // 375_E_N:375_E_CA:375_E_C:375_E_O [ 375_E_N:375_E_CA:375_E_C -- 375_E_CA:375_E_C:375_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.19433, 460, 462, 0, 0, 1,     // 375_E_N:375_E_CA:375_E_C:376_I_N [ 375_E_N:375_E_CA:375_E_C -- 375_E_CA:375_E_C:376_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.51234, 462, 463, 0, 0, 1,     // 375_E_CA:375_E_C:376_I_N:376_I_CA [ 375_E_CA:375_E_C:376_I_N -- 375_E_C:376_I_N:376_I_CA ] 
        [ [ 0.30327250802450434, 0.7454100770303694, 0.5936241259737559, 0.0 ], [ -0.33912433729307007, 0.6666061933867542, -0.6638003214777447, 0.0 ], [ -0.8905169676834142, 0.0, 0.45495003051756905, 15.376396220905745 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.88897, 463, 464, 0, 0, 1,     // 375_E_C:376_I_N:376_I_CA:376_I_C [ 375_E_C:376_I_N:376_I_CA -- 376_I_N:376_I_CA:376_I_C ] 
        [ [ -0.650583296862236, -0.7530323041160224, 0.09840590836722142, 7.9375535807794355 ], [ 0.7493582259974184, -0.6575772535632712, -0.07781005543107719, -8.875903771642113 ], [ 0.12330297228317588, 0.023119354530342022, 0.9920997290959386, 21.45969049046866 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 375, ' ') E sidechain

     ],
     [ // 116 : (' ', 376, ' ') I backbone
      [ 136.96491, 464, 465, 0, 0, 1,     // 376_I_N:376_I_CA:376_I_C:376_I_O [ 376_I_N:376_I_CA:376_I_C -- 376_I_CA:376_I_C:376_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.23116, 464, 466, 0, 0, 1,     // 376_I_N:376_I_CA:376_I_C:377_F_N [ 376_I_N:376_I_CA:376_I_C -- 376_I_CA:376_I_C:377_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.83320, 466, 467, 0, 0, 1,     // 376_I_CA:376_I_C:377_F_N:377_F_CA [ 376_I_CA:376_I_C:377_F_N -- 376_I_C:377_F_N:377_F_CA ] 
        [ [ 0.3505464341578529, 0.6590985862190877, 0.6653617445744938, 0.0 ], [ -0.3072171358810305, 0.7520565494987727, -0.58311969421128, 0.0 ], [ -0.8847230238443389, 0.0, 0.4661171216333181, 15.34578640896285 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.51179, 467, 468, 0, 0, 1,     // 376_I_C:377_F_N:377_F_CA:377_F_C [ 376_I_C:377_F_N:377_F_CA -- 377_F_N:377_F_CA:377_F_C ] 
        [ [ -0.7340274856977412, -0.6718810710432582, 0.09889123628520675, 8.887292142109693 ], [ 0.6712159513912076, -0.7399033056386028, -0.0448580528229261, -7.7887782375998915 ], [ 0.10330922920169527, 0.03345033152051886, 0.994086655419395, 21.57175155675439 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 376, ' ') I sidechain

     ],
     [ // 117 : (' ', 377, ' ') F backbone
      [ 129.80756, 468, 469, 0, 0, 1,     // 377_F_N:377_F_CA:377_F_C:377_F_O [ 377_F_N:377_F_CA:377_F_C -- 377_F_CA:377_F_C:377_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.20240, 468, 470, 0, 0, 1,     // 377_F_N:377_F_CA:377_F_C:378_D_N [ 377_F_N:377_F_CA:377_F_C -- 377_F_CA:377_F_C:378_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.04594, 470, 471, 0, 0, 1,     // 377_F_CA:377_F_C:378_D_N:378_D_CA [ 377_F_CA:377_F_C:378_D_N -- 377_F_C:378_D_N:378_D_CA ] 
        [ [ 0.3203400601143959, 0.7217891809609377, 0.6135166046111905, 0.0 ], [ -0.3340754992605416, 0.6921129808403673, -0.6398227743258919, 0.0 ], [ -0.8864399622533525, 0.0, 0.4628435948787397, 15.37247193710422 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.26148, 471, 472, 0, 0, 1,     // 377_F_C:378_D_N:378_D_CA:378_D_C [ 377_F_C:378_D_N:378_D_CA -- 378_D_N:378_D_CA:378_D_C ] 
        [ [ -0.6743800808257829, -0.7322924473850938, 0.09465346315988625, 8.191334628699657 ], [ 0.730638865201402, -0.6803191635133808, -0.05772940683621514, -8.542560067934579 ], [ 0.1066693734987602, 0.03022593686229073, 0.9938350152304852, 21.552103813253286 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 377, ' ') F sidechain

     ],
     [ // 118 : (' ', 378, ' ') D backbone
      [ 136.45095, 472, 473, 0, 0, 1,     // 378_D_N:378_D_CA:378_D_C:378_D_O [ 378_D_N:378_D_CA:378_D_C -- 378_D_CA:378_D_C:378_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.24514, 472, 474, 0, 0, 1,     // 378_D_N:378_D_CA:378_D_C:379_M_N [ 378_D_N:378_D_CA:378_D_C -- 378_D_CA:378_D_C:379_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.75631, 474, 475, 0, 0, 1,     // 378_D_CA:378_D_C:379_M_N:379_M_CA [ 378_D_CA:378_D_C:379_M_N -- 378_D_C:379_M_N:379_M_CA ] 
        [ [ 0.37982796071721947, 0.6190273269544536, 0.6874124589655197, 0.0 ], [ -0.2993799772640483, 0.7853694471289445, -0.5418177375554938, 0.0 ], [ -0.8752727286227855, 0.0, 0.48362966258204604, 15.425204627648204 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.85080, 475, 476, 0, 0, 1,     // 378_D_C:379_M_N:379_M_CA:379_M_C [ 378_D_C:379_M_N:379_M_CA -- 379_M_N:379_M_CA:379_M_C ] 
        [ [ -0.7379765177754717, -0.6619791531325949, 0.13105060102815816, 9.217473199420523 ], [ 0.6704226249421651, -0.7413525783645328, 0.030493581585576884, -7.265202150109608 ], [ 0.07696858565444688, 0.11036283509457248, 0.9909065957254711, 21.910166305537146 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 378, ' ') D sidechain

     ],
     [ // 119 : (' ', 379, ' ') M backbone
      [ 136.98202, 476, 477, 0, 0, 1,     // 379_M_N:379_M_CA:379_M_C:379_M_O [ 379_M_N:379_M_CA:379_M_C -- 379_M_CA:379_M_C:379_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.26482, 476, 478, 0, 0, 1,     // 379_M_N:379_M_CA:379_M_C:380_L_N [ 379_M_N:379_M_CA:379_M_C -- 379_M_CA:379_M_C:380_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.32596, 478, 479, 0, 0, 1,     // 379_M_CA:379_M_C:380_L_N:380_L_CA [ 379_M_CA:379_M_C:380_L_N -- 379_M_C:380_L_N:380_L_CA ] 
        [ [ 0.34338464172777067, 0.6725582252988921, 0.6555550483432305, 0.0 ], [ -0.31207076812436063, 0.7400442105596157, -0.5957737843338416, 0.0 ], [ -0.8858322773007101, 0.0, 0.46400557808310655, 15.39645879174874 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.78946, 479, 480, 0, 0, 1,     // 379_M_C:380_L_N:380_L_CA:380_L_C [ 379_M_C:380_L_N:380_L_CA -- 380_L_N:380_L_CA:380_L_C ] 
        [ [ -0.7175284935039716, -0.6878460983854315, 0.1096385240043851, 8.768563600833511 ], [ 0.6874511751610747, -0.7246790910098015, -0.04744572502202027, -7.96894224648755 ], [ 0.11208810275661879, 0.041327472571486334, 0.9928385051111156, 21.602897797665754 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 379, ' ') M sidechain

     ],
     [ // 120 : (' ', 380, ' ') L backbone
      [ 133.01308, 480, 481, 0, 0, 1,     // 380_L_N:380_L_CA:380_L_C:380_L_O [ 380_L_N:380_L_CA:380_L_C -- 380_L_CA:380_L_C:380_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.59280, 480, 482, 0, 0, 1,     // 380_L_N:380_L_CA:380_L_C:381_L_N [ 380_L_N:380_L_CA:380_L_C -- 380_L_CA:380_L_C:381_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.86256, 482, 483, 0, 0, 1,     // 380_L_CA:380_L_C:381_L_N:381_L_CA [ 380_L_CA:380_L_C:381_L_N -- 380_L_C:381_L_N:381_L_CA ] 
        [ [ 0.3478294494804224, 0.6895284893264074, 0.6352677675448273, 0.0 ], [ -0.33115012777949554, 0.7242585604652819, -0.6048050350851762, 0.0 ], [ -0.8771284210112964, 0.0, 0.48025590371616456, 15.329451569105467 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.37221, 483, 484, 0, 0, 1,     // 380_L_C:381_L_N:381_L_CA:381_L_C [ 380_L_C:381_L_N:381_L_CA -- 381_L_N:381_L_CA:381_L_C ] 
        [ [ -0.687826584723471, -0.7179049439805475, 0.10727106206210793, 8.501383665724905 ], [ 0.7221869711697635, -0.6916961368203121, 0.001559801429412631, -8.093720331652866 ], [ 0.07307919006316874, 0.07854263629483714, 0.9942285885366444, 21.756410253000052 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 380, ' ') L sidechain

     ],
     [ // 121 : (' ', 381, ' ') L backbone
      [ 134.27995, 484, 485, 0, 0, 1,     // 381_L_N:381_L_CA:381_L_C:381_L_O [ 381_L_N:381_L_CA:381_L_C -- 381_L_CA:381_L_C:381_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.91287, 484, 486, 0, 0, 1,     // 381_L_N:381_L_CA:381_L_C:382_A_N [ 381_L_N:381_L_CA:381_L_C -- 381_L_CA:381_L_C:382_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.52673, 486, 487, 0, 0, 1,     // 381_L_CA:381_L_C:382_A_N:382_A_CA [ 381_L_CA:381_L_C:382_A_N -- 381_L_C:382_A_N:382_A_CA ] 
        [ [ 0.3425489793192806, 0.693563726001281, 0.6337424987662883, 0.0 ], [ -0.32979053773274963, 0.7203952789789921, -0.6101383804078578, 0.0 ], [ -0.8797149526916448, 0.0, 0.4755014216705741, 15.354053269568158 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.30317, 487, 488, 0, 0, 1,     // 381_L_C:382_A_N:382_A_CA:382_A_C [ 381_L_C:382_A_N:382_A_CA -- 382_A_N:382_A_CA:382_A_C ] 
        [ [ -0.700784088680867, -0.7076997462017176, 0.0897927072637375, 8.471393360745145 ], [ 0.7077758068197988, -0.7054926745410957, -0.036517029581203404, -8.155871248945614 ], [ 0.0891911897684381, 0.037962552533729925, 0.9952907998533967, 21.710198856153607 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 381, ' ') L sidechain

     ],
     [ // 122 : (' ', 382, ' ') A backbone
      [ 129.76411, 488, 489, 0, 0, 1,     // 382_A_N:382_A_CA:382_A_C:382_A_O [ 382_A_N:382_A_CA:382_A_C -- 382_A_CA:382_A_C:382_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.49224, 488, 490, 0, 0, 1,     // 382_A_N:382_A_CA:382_A_C:383_T_N [ 382_A_N:382_A_CA:382_A_C -- 382_A_CA:382_A_C:383_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.38864, 490, 491, 0, 0, 1,     // 382_A_CA:382_A_C:383_T_N:383_T_CA [ 382_A_CA:382_A_C:383_T_N -- 382_A_C:383_T_N:383_T_CA ] 
        [ [ 0.32463384442567345, 0.7131554895136409, 0.6213067799645894, 0.0 ], [ -0.3302602929608567, 0.7010058828401935, -0.6320750676279281, 0.0 ], [ -0.8863075120672377, 0.0, 0.4630971756048435, 15.388441535115268 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.87444, 491, 492, 0, 0, 1,     // 382_A_C:383_T_N:383_T_CA:383_T_C [ 382_A_C:383_T_N:383_T_CA -- 383_T_N:383_T_CA:383_T_C ] 
        [ [ -0.6872518097490691, -0.7220021420425013, 0.07998660439516789, 8.297023254114547 ], [ 0.7208600975045377, -0.6914417965203805, -0.047633621008873896, -8.440824571001526 ], [ 0.0896976578422075, 0.02492285920010775, 0.9956571605060217, 21.572710174718946 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 382, ' ') A sidechain

     ],
     [ // 123 : (' ', 383, ' ') T backbone
      [ 136.35335, 492, 493, 0, 0, 1,     // 383_T_N:383_T_CA:383_T_C:383_T_O [ 383_T_N:383_T_CA:383_T_C -- 383_T_CA:383_T_C:383_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -36.71769, 492, 494, 0, 0, 1,     // 383_T_N:383_T_CA:383_T_C:384_T_N [ 383_T_N:383_T_CA:383_T_C -- 383_T_CA:383_T_C:384_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 171.59504, 494, 495, 0, 0, 1,     // 383_T_CA:383_T_C:384_T_N:384_T_CA [ 383_T_CA:383_T_C:384_T_N -- 383_T_C:384_T_N:384_T_CA ] 
        [ [ 0.38337177917520115, 0.5978726983292789, 0.703970393926128, 0.0 ], [ -0.28594071202846255, 0.8015910656890252, -0.525061589351372, 0.0 ], [ -0.8782163674953324, 0.0, 0.47826353808679933, 15.370988304936976 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.30390, 495, 496, 0, 0, 1,     // 383_T_C:384_T_N:384_T_CA:384_T_C [ 383_T_C:384_T_N:384_T_CA -- 384_T_N:384_T_CA:384_T_C ] 
        [ [ -0.747371822251355, -0.6474882894795178, 0.14897742879905393, 9.419497490114551 ], [ 0.6588513012971786, -0.7511861627420829, 0.040426620979736726, -7.025602732335645 ], [ 0.08573401940713668, 0.12836769021725516, 0.988013873396818, 21.770408274144728 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 383, ' ') T sidechain

     ],
     [ // 124 : (' ', 384, ' ') T backbone
      [ 141.00506, 496, 497, 0, 0, 1,     // 384_T_N:384_T_CA:384_T_C:384_T_O [ 384_T_N:384_T_CA:384_T_C -- 384_T_CA:384_T_C:384_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.55233, 496, 498, 0, 0, 1,     // 384_T_N:384_T_CA:384_T_C:385_S_N [ 384_T_N:384_T_CA:384_T_C -- 384_T_CA:384_T_C:385_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 170.19438, 498, 499, 0, 0, 1,     // 384_T_CA:384_T_C:385_S_N:385_S_CA [ 384_T_CA:384_T_C:385_S_N -- 384_T_C:385_S_N:385_S_CA ] 
        [ [ 0.3919596409967604, 0.581446240645449, 0.7129431317215738, 0.0 ], [ -0.2801225544463272, 0.8135848260871601, -0.5095204463524402, 0.0 ], [ -0.8762984618953494, 0.0, 0.4817686225563521, 15.36228127491384 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.35026, 499, 500, 0, 0, 1,     // 384_T_C:385_S_N:385_S_CA:385_S_C [ 384_T_C:385_S_N:385_S_CA -- 385_S_N:385_S_CA:385_S_C ] 
        [ [ -0.7561649397222548, -0.6397051343867002, 0.13781119321061164, 9.555789019599574 ], [ 0.6521748630579227, -0.7539927499028446, 0.07850401957432655, -6.8292542138095325 ], [ 0.05368921614457059, 0.14923898328935492, 0.9873424906973979, 21.819569533368487 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 384, ' ') T sidechain

     ],
     [ // 125 : (' ', 385, ' ') S backbone
      [ 131.57113, 500, 501, 0, 0, 1,     // 385_S_N:385_S_CA:385_S_C:385_S_O [ 385_S_N:385_S_CA:385_S_C -- 385_S_CA:385_S_C:385_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.77361, 500, 502, 0, 0, 1,     // 385_S_N:385_S_CA:385_S_C:386_R_N [ 385_S_N:385_S_CA:385_S_C -- 385_S_CA:385_S_C:386_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.32836, 502, 503, 0, 0, 1,     // 385_S_CA:385_S_C:386_R_N:386_R_CA [ 385_S_CA:385_S_C:386_R_N -- 385_S_C:386_R_N:386_R_CA ] 
        [ [ 0.3454056491763545, 0.6661890016701485, 0.6609743955485659, 0.0 ], [ -0.3085421079275227, 0.7457829537162477, -0.5904316671570979, 0.0 ], [ -0.8862825199408494, 0.0, 0.46314500412645926, 15.214190862467328 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.32715, 503, 504, 0, 0, 1,     // 385_S_C:386_R_N:386_R_CA:386_R_C [ 385_S_C:386_R_N:386_R_CA -- 386_R_N:386_R_CA:386_R_C ] 
        [ [ -0.7147356919855002, -0.6921074244574557, 0.10069857800820836, 8.84131694206817 ], [ 0.6958551601166092, -0.7181760690320121, 0.002954658832179726, -7.897724234290774 ], [ 0.07027436759655813, 0.07218342524841574, 0.9949125923307586, 21.409304923155062 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 385, ' ') S sidechain

     ],
     [ // 126 : (' ', 386, ' ') R backbone
      [ 142.34992, 504, 505, 0, 0, 1,     // 386_R_N:386_R_CA:386_R_C:386_R_O [ 386_R_N:386_R_CA:386_R_C -- 386_R_CA:386_R_C:386_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -34.45649, 504, 506, 0, 0, 1,     // 386_R_N:386_R_CA:386_R_C:387_F_N [ 386_R_N:386_R_CA:386_R_C -- 386_R_CA:386_R_C:387_F_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.10233, 506, 507, 0, 0, 1,     // 386_R_CA:386_R_C:387_F_N:387_F_CA [ 386_R_CA:386_R_C:387_F_N -- 386_R_C:387_F_N:387_F_CA ] 
        [ [ 0.40577914128831655, 0.5657802073993148, 0.7177994465103082, 0.0 ], [ -0.27843079152879174, 0.8245560968880095, -0.49252770217873243, 0.0 ], [ -0.8705283354515044, 0.0, 0.4921182959167775, 15.327325757874437 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.99272, 507, 508, 0, 0, 1,     // 386_R_C:387_F_N:387_F_CA:387_F_C [ 386_R_C:387_F_N:387_F_CA -- 387_F_N:387_F_CA:387_F_C ] 
        [ [ -0.7905922406774482, -0.6044801820155167, 0.09781420413774082, 9.595705240117834 ], [ 0.609899454795691, -0.7915823759814866, 0.037682848557708384, -6.584221645302963 ], [ 0.054649464961064334, 0.08944859745123018, 0.9944910177540455, 21.906074367732366 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 386, ' ') R sidechain

     ],
     [ // 127 : (' ', 387, ' ') F backbone
      [ 134.14440, 508, 509, 0, 0, 1,     // 387_F_N:387_F_CA:387_F_C:387_F_O [ 387_F_N:387_F_CA:387_F_C -- 387_F_CA:387_F_C:387_F_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.57768, 508, 510, 0, 0, 1,     // 387_F_N:387_F_CA:387_F_C:388_R_N [ 387_F_N:387_F_CA:387_F_C -- 387_F_CA:387_F_C:388_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.13693, 510, 511, 0, 0, 1,     // 387_F_CA:387_F_C:388_R_N:388_R_CA [ 387_F_CA:387_F_C:388_R_N -- 387_F_C:388_R_N:388_R_CA ] 
        [ [ 0.3480162449771242, 0.6765891563982366, 0.6489312803959636, 0.0 ], [ -0.31976721523817847, 0.7363607223666417, -0.5962565006895085, 0.0 ], [ -0.8812681891971608, 0.0, 0.4726165239484939, 15.407290668592633 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.72234, 511, 512, 0, 0, 1,     // 387_F_C:388_R_N:388_R_CA:388_R_C [ 387_F_C:388_R_N:388_R_CA -- 388_R_N:388_R_CA:388_R_C ] 
        [ [ -0.6931617357188652, -0.7086001854233202, 0.13195675561806655, 8.732731123101134 ], [ 0.7142696218697014, -0.6998440320397026, -0.006102302231560138, -8.023881508293881 ], [ 0.09667324039941803, 0.09002281953175967, 0.9912366904802448, 21.76733702890761 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 387, ' ') F sidechain

     ],
     [ // 128 : (' ', 388, ' ') R backbone
      [ 128.96686, 512, 513, 0, 0, 1,     // 388_R_N:388_R_CA:388_R_C:388_R_O [ 388_R_N:388_R_CA:388_R_C -- 388_R_CA:388_R_C:388_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.89982, 512, 514, 0, 0, 1,     // 388_R_N:388_R_CA:388_R_C:389_E_N [ 388_R_N:388_R_CA:388_R_C -- 388_R_CA:388_R_C:389_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.13339, 514, 515, 0, 0, 1,     // 388_R_CA:388_R_C:389_E_N:389_E_CA [ 388_R_CA:388_R_C:389_E_N -- 388_R_C:389_E_N:389_E_CA ] 
        [ [ 0.30269583766294117, 0.7649193985463456, 0.5685714938237998, 0.0 ], [ -0.3594605951921811, 0.6441260076471814, -0.6751960950543079, 0.0 ], [ -0.8827022773085007, 0.0, 0.4699326437207642, 15.28014578075692 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -56.84536, 515, 516, 0, 0, 1,     // 388_R_C:389_E_N:389_E_CA:389_E_C [ 388_R_C:389_E_N:389_E_CA -- 389_E_N:389_E_CA:389_E_C ] 
        [ [ -0.626406519669156, -0.7743731594541474, 0.08922489581376396, 7.579969310241528 ], [ 0.7722511321106317, -0.6320755707059296, -0.06409884453421213, -9.00144614090103 ], [ 0.10603330170198003, 0.028752092685142262, 0.9939468074783488, 21.545101122306185 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 388, ' ') R sidechain

     ],
     [ // 129 : (' ', 389, ' ') E backbone
      [ 134.98825, 516, 517, 0, 0, 1,     // 389_E_N:389_E_CA:389_E_C:389_E_O [ 389_E_N:389_E_CA:389_E_C -- 389_E_CA:389_E_C:389_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.42988, 516, 518, 0, 0, 1,     // 389_E_N:389_E_CA:389_E_C:390_L_N [ 389_E_N:389_E_CA:389_E_C -- 389_E_CA:389_E_C:390_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.03714, 518, 519, 0, 0, 1,     // 389_E_CA:389_E_C:390_L_N:390_L_CA [ 389_E_CA:389_E_C:390_L_N -- 389_E_C:390_L_N:390_L_CA ] 
        [ [ 0.3198415498204168, 0.687466300739426, 0.6519903897728271, 0.0 ], [ -0.30277521054140355, 0.7262162800073054, -0.6172010098290182, 0.0 ], [ -0.8977909305011288, 0.0, 0.4404218944488539, 15.304332028868759 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -73.08800, 519, 520, 0, 0, 1,     // 389_E_C:390_L_N:390_L_CA:390_L_C [ 389_E_C:390_L_N:390_L_CA -- 390_L_N:390_L_CA:390_L_C ] 
        [ [ -0.6926773596072175, -0.7079267849494447, 0.13797732653851547, 8.70356570833206 ], [ 0.7087642274681687, -0.7035552191750165, -0.05160739707681811, -8.239154485340611 ], [ 0.1336089269061893, 0.06204611770882258, 0.9890899523947453, 21.1836225181953 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 389, ' ') E sidechain

     ],
     [ // 130 : (' ', 390, ' ') L backbone
      [ 171.26038, 520, 521, 0, 0, 1,     // 390_L_N:390_L_CA:390_L_C:390_L_O [ 390_L_N:390_L_CA:390_L_C -- 390_L_CA:390_L_C:390_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -9.44486, 520, 522, 0, 0, 1,     // 390_L_N:390_L_CA:390_L_C:391_K_N [ 390_L_N:390_L_CA:390_L_C -- 390_L_CA:390_L_C:391_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -171.05335, 522, 523, 0, 0, 1,     // 390_L_CA:390_L_C:391_K_N:391_K_CA [ 390_L_CA:390_L_C:391_K_N -- 390_L_C:391_K_N:391_K_CA ] 
        [ [ 0.501546405951398, 0.16409833563910034, 0.8494250637447147, 0.0 ], [ -0.0834339624797722, 0.9864439853536932, -0.1413047687250745, 0.0 ], [ -0.8610981225053039, 0.0, 0.508438809905224, 15.340933671363157 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  46.18805, 523, 524, 0, 0, 1,     // 390_L_C:391_K_N:391_K_CA:391_K_C [ 390_L_C:391_K_N:391_K_CA -- 391_K_N:391_K_CA:391_K_C ] 
        [ [ -0.994971362753535, -0.08410395510174762, 0.05439220566052376, 11.37547838423799 ], [ 0.07588101818991613, -0.9874176922236836, -0.13873850281056369, -1.8923497914407226 ], [ 0.06537628299953367, -0.13391350126069582, 0.9888342205856685, 22.14993300355475 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 390, ' ') L sidechain

     ],
     [ // 131 : (' ', 391, ' ') K backbone
      [ -118.65782, 524, 525, 0, 0, 1,     // 391_K_N:391_K_CA:391_K_C:391_K_O [ 391_K_N:391_K_CA:391_K_C -- 391_K_CA:391_K_C:391_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  54.21606, 524, 526, 0, 0, 1,     // 391_K_N:391_K_CA:391_K_C:392_L_N [ 391_K_N:391_K_CA:391_K_C -- 391_K_CA:391_K_C:392_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -175.38058, 526, 527, 0, 0, 1,     // 391_K_CA:391_K_C:392_L_N:392_L_CA [ 391_K_CA:391_K_C:392_L_N -- 391_K_C:392_L_N:392_L_CA ] 
        [ [ 0.2590188270721803, -0.8112277003202923, 0.5242316906246721, 0.0 ], [ 0.35935065856705484, 0.5847303808021693, 0.7272946349000342, 0.0 ], [ -0.8965357502127709, 0.0, 0.4429713857467817, 15.339036995602894 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -56.39511, 527, 528, 0, 0, 1,     // 391_K_C:392_L_N:392_L_CA:392_L_C [ 391_K_C:392_L_N:392_L_CA -- 392_L_N:392_L_CA:392_L_C ] 
        [ [ -0.5433065181472845, 0.8294530917823226, 0.12971351460592068, 7.005113435405221 ], [ -0.8300275258093437, -0.5538899971479196, 0.06527003491876468, 9.718568162801157 ], [ 0.12598545049922957, -0.0722041521813208, 0.9894009433340373, 21.25829914874982 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 391, ' ') K sidechain

     ],
     [ // 132 : (' ', 392, ' ') L backbone
      [ -59.51227, 528, 529, 0, 0, 1,     // 392_L_N:392_L_CA:392_L_C:392_L_O [ 392_L_N:392_L_CA:392_L_C -- 392_L_CA:392_L_C:392_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 120.27107, 528, 530, 0, 0, 1,     // 392_L_N:392_L_CA:392_L_C:393_Q_N [ 392_L_N:392_L_CA:392_L_C -- 392_L_CA:392_L_C:393_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.02130, 530, 531, 0, 0, 1,     // 392_L_CA:392_L_C:393_Q_N:393_Q_CA [ 392_L_CA:392_L_C:393_Q_N -- 392_L_C:393_Q_N:393_Q_CA ] 
        [ [ -0.23328832402865723, -0.8636502224121513, -0.4468611095176463, 0.0 ], [ 0.3996883342403684, -0.5040915525253737, 0.7655984209392793, 0.0 ], [ -0.8864681569825618, 0.0, 0.46278959220788474, 15.362159447716966 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -79.47028, 531, 532, 0, 0, 1,     // 392_L_C:393_Q_N:393_Q_CA:393_Q_C [ 392_L_C:393_Q_N:393_Q_CA -- 393_Q_N:393_Q_CA:393_Q_C ] 
        [ [ 0.48643425143830027, 0.8711902271997356, -0.06640261334716666, -5.953186881658006 ], [ -0.8699380059766078, 0.4899905527413136, 0.05583120974603132, 10.199478941171582 ], [ 0.08117625752092705, 0.03060794432716436, 0.9962296768110059, 21.52754925367644 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 392, ' ') L sidechain

     ],
     [ // 133 : (' ', 393, ' ') Q backbone
      [ -35.39352, 532, 533, 0, 0, 1,     // 393_Q_N:393_Q_CA:393_Q_C:393_Q_O [ 393_Q_N:393_Q_CA:393_Q_C -- 393_Q_CA:393_Q_C:393_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 150.16934, 532, 534, 0, 0, 1,     // 393_Q_N:393_Q_CA:393_Q_C:394_H_N [ 393_Q_N:393_Q_CA:393_Q_C -- 393_Q_CA:393_Q_C:394_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.17469, 534, 535, 0, 0, 1,     // 393_Q_CA:393_Q_C:394_H_N:394_H_CA [ 393_Q_CA:393_Q_C:394_H_N -- 393_Q_C:394_H_N:394_H_CA ] 
        [ [ -0.39534445347606895, -0.49743826133046987, -0.7721774014241983, 0.0 ], [ 0.22669694279229405, -0.86749938107443, 0.44277908708983554, 0.0 ], [ -0.8901186770506142, 0.0, 0.45572880177323033, 15.285924207881873 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.06231, 535, 536, 0, 0, 1,     // 393_Q_C:394_H_N:394_H_CA:394_H_C [ 393_Q_C:394_H_N:394_H_CA -- 394_H_N:394_H_CA:394_H_C ] 
        [ [ 0.8677734894775052, 0.4916921581548064, -0.07216642272624542, -10.313358914663914 ], [ -0.4887272549031031, 0.8706747034100308, 0.055418689597754096, 5.913847823883751 ], [ 0.09008241379376768, -0.01282117197933895, 0.9958517842902952, 21.372731063104816 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 393, ' ') Q sidechain

     ],
     [ // 134 : (' ', 394, ' ') H backbone
      [ 132.48972, 536, 537, 0, 0, 1,     // 394_H_N:394_H_CA:394_H_C:394_H_O [ 394_H_N:394_H_CA:394_H_C -- 394_H_CA:394_H_C:394_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.34026, 536, 538, 0, 0, 1,     // 394_H_N:394_H_CA:394_H_C:395_K_N [ 394_H_N:394_H_CA:394_H_C -- 394_H_CA:394_H_C:395_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.70108, 538, 539, 0, 0, 1,     // 394_H_CA:394_H_C:395_K_N:395_K_CA [ 394_H_CA:394_H_C:395_K_N -- 394_H_C:395_K_N:395_K_CA ] 
        [ [ 0.3487584429581048, 0.6735320679295658, 0.6517070675817175, 0.0 ], [ -0.3177940228269184, 0.7391580030484165, -0.5938454391379641, 0.0 ], [ -0.8816884412993756, 0.0, 0.47183205961345775, 15.39696182755672 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.99300, 539, 540, 0, 0, 1,     // 394_H_C:395_K_N:395_K_CA:395_K_C [ 394_H_C:395_K_N:395_K_CA -- 395_K_N:395_K_CA:395_K_C ] 
        [ [ -0.7152356479099072, -0.692485336262833, 0.0943505538926315, 8.735243839120612 ], [ 0.6940354004767174, -0.719645622133135, -0.0206165329222803, -7.959687675121505 ], [ 0.08217560978792772, 0.05073694517374953, 0.9953255404893507, 21.72122666502719 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 394, ' ') H sidechain

     ],
     [ // 135 : (' ', 395, ' ') K backbone
      [ 141.79850, 540, 541, 0, 0, 1,     // 395_K_N:395_K_CA:395_K_C:395_K_O [ 395_K_N:395_K_CA:395_K_C -- 395_K_CA:395_K_C:395_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.57864, 540, 542, 0, 0, 1,     // 395_K_N:395_K_CA:395_K_C:396_E_N [ 395_K_N:395_K_CA:395_K_C -- 395_K_CA:395_K_C:396_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.30331, 542, 543, 0, 0, 1,     // 395_K_CA:395_K_C:396_E_N:396_E_CA [ 395_K_CA:395_K_C:396_E_N -- 395_K_C:396_E_N:396_E_CA ] 
        [ [ 0.3889575958849897, 0.5818197330356681, 0.7142813079268371, 0.0 ], [ -0.27824696864010806, 0.813317771999545, -0.5109724319298321, 0.0 ], [ -0.878231525877977, 0.0, 0.4782357022996508, 15.369964984164415 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.15764, 543, 544, 0, 0, 1,     // 395_K_C:396_E_N:396_E_CA:396_E_C [ 395_K_C:396_E_N:396_E_CA -- 396_E_N:396_E_CA:396_E_C ] 
        [ [ -0.7840776378362472, -0.6117142763434916, 0.10501381796072738, 9.527821234034683 ], [ 0.6158573527531198, -0.7878036118744342, 0.009229852518595814, -6.815877628210403 ], [ 0.07708423253201536, 0.07191045289214913, 0.9944279299475635, 21.749166327776816 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 395, ' ') K sidechain

     ],
     [ // 136 : (' ', 396, ' ') E backbone
      [ 132.94130, 544, 545, 0, 0, 1,     // 396_E_N:396_E_CA:396_E_C:396_E_O [ 396_E_N:396_E_CA:396_E_C -- 396_E_CA:396_E_C:396_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.03346, 544, 546, 0, 0, 1,     // 396_E_N:396_E_CA:396_E_C:397_Y_N [ 396_E_N:396_E_CA:396_E_C -- 396_E_CA:396_E_C:397_Y_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.19298, 546, 547, 0, 0, 1,     // 396_E_CA:396_E_C:397_Y_N:397_Y_CA [ 396_E_CA:396_E_C:397_Y_N -- 396_E_C:397_Y_N:397_Y_CA ] 
        [ [ 0.3353063831499286, 0.6950783421760851, 0.6359526143170247, 0.0 ], [ -0.32418025089459573, 0.718934001308705, -0.6148504425404374, 0.0 ], [ -0.8845771839409098, 0.0, 0.4663938310603709, 15.357543769580431 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -53.97546, 547, 548, 0, 0, 1,     // 396_E_C:397_Y_N:397_Y_CA:397_Y_C [ 396_E_C:397_Y_N:397_Y_CA -- 397_Y_N:397_Y_CA:397_Y_C ] 
        [ [ -0.6741362558329481, -0.7254370974445471, 0.13885721523480424, 8.520716038558323 ], [ 0.7307216529756315, -0.6824447179092359, -0.01775029211997311, -8.237981744434414 ], [ 0.10763909347488412, 0.08949985837426942, 0.9901533219188283, 21.60645179690038 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 396, ' ') E sidechain

     ],
     [ // 137 : (' ', 397, ' ') Y backbone
      [ 128.85062, 548, 549, 0, 0, 1,     // 397_Y_N:397_Y_CA:397_Y_C:397_Y_O [ 397_Y_N:397_Y_CA:397_Y_C -- 397_Y_CA:397_Y_C:397_Y_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.02014, 548, 550, 0, 0, 1,     // 397_Y_N:397_Y_CA:397_Y_C:398_L_N [ 397_Y_N:397_Y_CA:397_Y_C -- 397_Y_CA:397_Y_C:398_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.96429, 550, 551, 0, 0, 1,     // 397_Y_CA:397_Y_C:398_L_N:398_L_CA [ 397_Y_CA:397_Y_C:398_L_N -- 397_Y_C:398_L_N:398_L_CA ] 
        [ [ 0.31423455049170834, 0.7433800420355205, 0.5904597872679745, 0.0 ], [ -0.34923968031027325, 0.6688692795332043, -0.6562358818241428, 0.0 ], [ -0.8827730699189074, 0.0, 0.4697996456213519, 15.362925356229043 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.25511, 551, 552, 0, 0, 1,     // 397_Y_C:398_L_N:398_L_CA:398_L_C [ 397_Y_C:398_L_N:398_L_CA -- 398_L_N:398_L_CA:398_L_C ] 
        [ [ -0.6655627699283151, -0.7435757507551857, 0.0641973689041119, 7.9103953373292395 ], [ 0.7402194450818723, -0.6686514788451583, -0.07057175752996292, -8.791598296284853 ], [ 0.09540111324315771, 0.0005502063655663984, 0.9954387599771879, 21.656835638243663 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 397, ' ') Y sidechain

     ],
     [ // 138 : (' ', 398, ' ') L backbone
      [ 132.34076, 552, 553, 0, 0, 1,     // 398_L_N:398_L_CA:398_L_C:398_L_O [ 398_L_N:398_L_CA:398_L_C -- 398_L_CA:398_L_C:398_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.80777, 552, 554, 0, 0, 1,     // 398_L_N:398_L_CA:398_L_C:399_C_N [ 398_L_N:398_L_CA:398_L_C -- 398_L_CA:398_L_C:399_C_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.96160, 554, 555, 0, 0, 1,     // 398_L_CA:398_L_C:399_C_N:399_C_CA [ 398_L_CA:398_L_C:399_C_N -- 398_L_C:399_C_N:399_C_CA ] 
        [ [ 0.32909949073431216, 0.7047304579857736, 0.6285288432407682, 0.0 ], [ -0.32689860443058255, 0.7094751451510911, -0.6243254926992615, 0.0 ], [ -0.8859067826920357, 0.0, 0.46386331217315097, 15.315902489612599 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.85154, 555, 556, 0, 0, 1,     // 398_L_C:399_C_N:399_C_CA:399_C_C [ 398_L_C:399_C_N:399_C_CA -- 399_C_N:399_C_CA:399_C_C ] 
        [ [ -0.6805100987053841, -0.7261574834349231, 0.09798629910067803, 8.39662638608685 ], [ 0.7287123878627583, -0.6846917590372081, -0.01324578732619025, -8.340473093447475 ], [ 0.07670893908368975, 0.06238993795471327, 0.9950996102434478, 21.51273334060144 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 398, ' ') L sidechain

     ],
     [ // 139 : (' ', 399, ' ') C backbone
      [ 128.54523, 556, 557, 0, 0, 1,     // 399_C_N:399_C_CA:399_C_C:399_C_O [ 399_C_N:399_C_CA:399_C_C -- 399_C_CA:399_C_C:399_C_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.11703, 556, 558, 0, 0, 1,     // 399_C_N:399_C_CA:399_C_C:400_V_N [ 399_C_N:399_C_CA:399_C_C -- 399_C_CA:399_C_C:400_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.31643, 558, 559, 0, 0, 1,     // 399_C_CA:399_C_C:400_V_N:400_V_CA [ 399_C_CA:399_C_C:400_V_N -- 399_C_C:400_V_N:400_V_CA ] 
        [ [ 0.3328459332007128, 0.7445099897770938, 0.5787214009122682, 0.0 ], [ -0.3711847219386294, 0.6676113203969145, -0.6453813036316056, 0.0 ], [ -0.8668537863740857, 0.0, 0.4985624464888135, 15.356518468812801 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.09382, 559, 560, 0, 0, 1,     // 399_C_C:400_V_N:400_V_CA:400_V_C [ 399_C_C:400_V_N:400_V_CA -- 400_V_N:400_V_CA:400_V_C ] 
        [ [ -0.6621662769779416, -0.7484279234969288, 0.03730234527810712, 7.751217465930577 ], [ 0.7480293636430669, -0.6631355003048226, -0.026521300179682898, -8.644039817792574 ], [ 0.044585791020459606, 0.010341739000181396, 0.9989520287148588, 22.034111412480364 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 399, ' ') C sidechain

     ],
     [ // 140 : (' ', 400, ' ') V backbone
      [ 136.54237, 560, 561, 0, 0, 1,     // 400_V_N:400_V_CA:400_V_C:400_V_O [ 400_V_N:400_V_CA:400_V_C -- 400_V_CA:400_V_C:400_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.79559, 560, 562, 0, 0, 1,     // 400_V_N:400_V_CA:400_V_C:401_K_N [ 400_V_N:400_V_CA:400_V_C -- 400_V_CA:400_V_C:401_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.77665, 562, 563, 0, 0, 1,     // 400_V_CA:400_V_C:401_K_N:401_K_CA [ 400_V_CA:400_V_C:401_K_N -- 400_V_C:401_K_N:401_K_CA ] 
        [ [ 0.3726318583720551, 0.6265437931320581, 0.6845351513354748, 0.0 ], [ -0.29955646236293304, 0.7793862170244562, -0.5502936039687468, 0.0 ], [ -0.8783003039865086, 0.0, 0.47810937662548164, 15.379954915964435 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.88290, 563, 564, 0, 0, 1,     // 400_V_C:401_K_N:401_K_CA:401_K_C [ 400_V_C:401_K_N:401_K_CA -- 401_K_N:401_K_CA:401_K_C ] 
        [ [ -0.7349658461331681, -0.6632464315467201, 0.1411714420777027, 9.166082805193481 ], [ 0.67004505126786, -0.7423199557634431, 0.0008441248785734967, -7.368557672028029 ], [ 0.10423451582468808, 0.09521162910013488, 0.9899848036181622, 21.781949214633975 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 400, ' ') V sidechain

     ],
     [ // 141 : (' ', 401, ' ') K backbone
      [ 134.29280, 564, 565, 0, 0, 1,     // 401_K_N:401_K_CA:401_K_C:401_K_O [ 401_K_N:401_K_CA:401_K_C -- 401_K_CA:401_K_C:401_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.95567, 564, 566, 0, 0, 1,     // 401_K_N:401_K_CA:401_K_C:402_A_N [ 401_K_N:401_K_CA:401_K_C -- 401_K_CA:401_K_C:402_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.07998, 566, 567, 0, 0, 1,     // 401_K_CA:401_K_C:402_A_N:402_A_CA [ 401_K_CA:401_K_C:402_A_N -- 401_K_C:402_A_N:402_A_CA ] 
        [ [ 0.3483184605836771, 0.681432346720359, 0.6436802054279763, 0.0 ], [ -0.3243087753940719, 0.7318811084070858, -0.5993111557106316, 0.0 ], [ -0.8794873894599142, 0.0, 0.4759221908894196, 15.371156843005616 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.19569, 567, 568, 0, 0, 1,     // 401_K_C:402_A_N:402_A_CA:402_A_C [ 401_K_C:402_A_N:402_A_CA -- 402_A_N:402_A_CA:402_A_C ] 
        [ [ -0.6902589988020704, -0.7137237672600169, 0.11891551043044332, 8.581100888682835 ], [ 0.719426249792428, -0.6945286167060322, 0.0074747365176719375, -7.98960329601355 ], [ 0.0772553278574787, 0.09071044385612369, 0.9928762408744901, 21.715823496953305 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 401, ' ') K sidechain

     ],
     [ // 142 : (' ', 402, ' ') A backbone
      [ 140.63770, 568, 569, 0, 0, 1,     // 402_A_N:402_A_CA:402_A_C:402_A_O [ 402_A_N:402_A_CA:402_A_C -- 402_A_CA:402_A_C:402_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.26351, 568, 570, 0, 0, 1,     // 402_A_N:402_A_CA:402_A_C:403_M_N [ 402_A_N:402_A_CA:402_A_C -- 402_A_CA:402_A_C:403_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.26678, 570, 571, 0, 0, 1,     // 402_A_CA:402_A_C:403_M_N:403_M_CA [ 402_A_CA:402_A_C:403_M_N -- 402_A_C:403_M_N:403_M_CA ] 
        [ [ 0.3927817842029261, 0.57733769499126, 0.7158237603911627, 0.0 ], [ -0.27772958991645547, 0.8165054720834263, -0.5061472996496756, 0.0 ], [ -0.8766919327125141, 0.0, 0.48105223740961495, 15.437047425019237 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.97233, 571, 572, 0, 0, 1,     // 402_A_C:403_M_N:403_M_CA:403_M_C [ 402_A_C:403_M_N:403_M_CA -- 403_M_N:403_M_CA:403_M_C ] 
        [ [ -0.7707233210828928, -0.6249399979878075, 0.12423993421583691, 9.597023920775667 ], [ 0.6339389711081288, -0.7717081623410796, 0.05087133855641064, -6.785899003295946 ], [ 0.0640854371079998, 0.11796826306745828, 0.9909472971149997, 21.886497813533786 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 402, ' ') A sidechain

     ],
     [ // 143 : (' ', 403, ' ') M backbone
      [ 136.98250, 572, 573, 0, 0, 1,     // 403_M_N:403_M_CA:403_M_C:403_M_O [ 403_M_N:403_M_CA:403_M_C -- 403_M_CA:403_M_C:403_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.16003, 572, 574, 0, 0, 1,     // 403_M_N:403_M_CA:403_M_C:404_I_N [ 403_M_N:403_M_CA:403_M_C -- 403_M_CA:403_M_C:404_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.85724, 574, 575, 0, 0, 1,     // 403_M_CA:403_M_C:404_I_N:404_I_CA [ 403_M_CA:403_M_C:404_I_N -- 403_M_C:404_I_N:404_I_CA ] 
        [ [ 0.3826346592613615, 0.6449247540503186, 0.6615608657901991, 0.0 ], [ -0.32289411998686085, 0.7642460739926219, -0.5582717310277985, 0.0 ], [ -0.8656385532136681, 0.0, 0.5006694470308204, 15.296741202740508 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -55.50761, 575, 576, 0, 0, 1,     // 403_M_C:404_I_N:404_I_CA:404_I_C [ 403_M_C:404_I_N:404_I_CA -- 404_I_N:404_I_CA:404_I_C ] 
        [ [ -0.7371591323096215, -0.6708818927946311, 0.08070873294024225, 8.863914585816103 ], [ 0.6736991830106441, -0.7389227421178509, 0.01107212679780314, -7.479996468042922 ], [ 0.052209428873688804, 0.06253532682675614, 0.9966762304958185, 22.00495422393871 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 403, ' ') M sidechain

     ],
     [ // 144 : (' ', 404, ' ') I backbone
      [ 129.47525, 576, 577, 0, 0, 1,     // 404_I_N:404_I_CA:404_I_C:404_I_O [ 404_I_N:404_I_CA:404_I_C -- 404_I_CA:404_I_C:404_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.54646, 576, 578, 0, 0, 1,     // 404_I_N:404_I_CA:404_I_C:405_L_N [ 404_I_N:404_I_CA:404_I_C -- 404_I_CA:404_I_C:405_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.74765, 578, 579, 0, 0, 1,     // 404_I_CA:404_I_C:405_L_N:405_L_CA [ 404_I_CA:404_I_C:405_L_N -- 404_I_C:405_L_N:405_L_CA ] 
        [ [ 0.3270809993902617, 0.7259323213415708, 0.6050126318263991, 0.0 ], [ -0.34523169282125044, 0.6877661410913148, -0.6385865747415423, 0.0 ], [ -0.8796778376824332, 0.0, 0.4755700809453414, 15.394044992325451 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.29086, 579, 580, 0, 0, 1,     // 404_I_C:405_L_N:405_L_CA:405_L_C [ 404_I_C:405_L_N:405_L_CA -- 405_L_N:405_L_CA:405_L_C ] 
        [ [ -0.669937408772595, -0.7382260265277044, 0.07877945216982703, 8.109119260793696 ], [ 0.7381962887117776, -0.6736669118805039, -0.03520129499810196, -8.559118306817078 ], [ 0.0790576223979636, 0.034572034862045166, 0.9962703783342545, 21.768216869711388 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 404, ' ') I sidechain

     ],
     [ // 145 : (' ', 405, ' ') L backbone
      [ 141.29867, 580, 581, 0, 0, 1,     // 405_L_N:405_L_CA:405_L_C:405_L_O [ 405_L_N:405_L_CA:405_L_C -- 405_L_CA:405_L_C:405_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.18840, 580, 582, 0, 0, 1,     // 405_L_N:405_L_CA:405_L_C:406_L_N [ 405_L_N:405_L_CA:405_L_C -- 405_L_CA:405_L_C:406_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -173.13457, 582, 583, 0, 0, 1,     // 405_L_CA:405_L_C:406_L_N:406_L_CA [ 405_L_CA:405_L_C:406_L_N -- 405_L_C:406_L_N:406_L_CA ] 
        [ [ 0.35068482902035175, 0.6453030080040858, 0.6786782584965025, 0.0 ], [ -0.296229952694139, 0.7639267162895135, -0.5732920610526099, 0.0 ], [ -0.8884075448924299, 0.0, 0.4590555894205027, 15.377389340042024 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -89.64834, 583, 584, 0, 0, 1,     // 405_L_C:406_L_N:406_L_CA:406_L_C [ 405_L_C:406_L_N:406_L_CA -- 406_L_N:406_L_CA:406_L_C ] 
        [ [ -0.8005042253742721, -0.5987558848234767, 0.026160572378906854, 9.057308549654815 ], [ 0.5888216662684209, -0.7938597627659759, -0.15190695307993243, -7.650875832569186 ], [ 0.11172300788475445, -0.10619824598553365, 0.9880485322385633, 21.50372101896758 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 405, ' ') L sidechain

     ],
     [ // 146 : (' ', 406, ' ') L backbone
      [ 164.68903, 584, 585, 0, 0, 1,     // 406_L_N:406_L_CA:406_L_C:406_L_O [ 406_L_N:406_L_CA:406_L_C -- 406_L_CA:406_L_C:406_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -13.84335, 584, 586, 0, 0, 1,     // 406_L_N:406_L_CA:406_L_C:407_N_N [ 406_L_N:406_L_CA:406_L_C -- 406_L_CA:406_L_C:407_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -175.03136, 586, 587, 0, 0, 1,     // 406_L_CA:406_L_C:407_N_N:407_N_CA [ 406_L_CA:406_L_C:407_N_N -- 406_L_C:407_N_N:407_N_CA ] 
        [ [ 0.4670762585030024, 0.2392680904885779, 0.8512288467954937, 0.0 ], [ -0.11509968257119452, 0.9709535420780696, -0.20976482593147153, 0.0 ], [ -0.8766936932674072, 0.0, 0.4810490288787134, 15.42482292785974 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -95.32556, 587, 588, 0, 0, 1,     // 406_L_C:407_N_N:407_N_CA:407_N_C [ 406_L_C:407_N_N:407_N_CA -- 407_N_N:407_N_CA:407_N_C ] 
        [ [ -0.9752446977268254, -0.19791528233053496, 0.09862717970081891, 11.355960935250776 ], [ 0.18852507695397033, -0.9772737821457178, -0.09692393971617047, -2.7984027771556015 ], [ 0.11556848582209828, -0.07593086164813977, 0.9903930680967827, 21.842338030234338 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 406, ' ') L sidechain

     ],
     [ // 147 : (' ', 407, ' ') N backbone
      [ -144.71781, 588, 589, 0, 0, 1,     // 407_N_N:407_N_CA:407_N_C:407_N_O [ 407_N_N:407_N_CA:407_N_C -- 407_N_CA:407_N_C:407_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  35.68152, 588, 590, 0, 0, 1,     // 407_N_N:407_N_CA:407_N_C:408_S_N [ 407_N_N:407_N_CA:407_N_C -- 407_N_CA:407_N_C:408_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 168.20673, 590, 591, 0, 0, 1,     // 407_N_CA:407_N_C:408_S_N:408_S_CA [ 407_N_CA:407_N_C:408_S_N -- 407_N_C:408_S_N:408_S_CA ] 
        [ [ 0.3772100496211699, -0.5832791924078338, 0.7193733121049597, 0.0 ], [ 0.2708684318222141, 0.8122717425246709, 0.5165703330023022, 0.0 ], [ -0.8856313404046681, 0.0, 0.46438898446564264, 15.460344279109044 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -71.32909, 591, 592, 0, 0, 1,     // 407_N_C:408_S_N:408_S_CA:408_S_C [ 407_N_C:408_S_N:408_S_CA -- 408_S_N:408_S_CA:408_S_C ] 
        [ [ -0.8685818288597863, 0.4938724300157565, -0.04068942669687905, 9.615129571426122 ], [ -0.4922653615618604, -0.8504861978539243, 0.18533224507987378, 6.904468932880266 ], [ 0.05692469043358835, 0.18100621572285086, 0.9818331474791085, 21.667358331541198 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 407, ' ') N sidechain

     ],
     [ // 148 : (' ', 408, ' ') S backbone
      [ -56.13708, 592, 593, 0, 0, 1,     // 408_S_N:408_S_CA:408_S_C:408_S_O [ 408_S_N:408_S_CA:408_S_C -- 408_S_CA:408_S_C:408_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 130.50792, 592, 594, 0, 0, 1,     // 408_S_N:408_S_CA:408_S_C:409_S_N [ 408_S_N:408_S_CA:408_S_C -- 408_S_CA:408_S_C:409_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 168.32051, 594, 595, 0, 0, 1,     // 408_S_CA:408_S_C:409_S_N:409_S_CA [ 408_S_CA:408_S_C:409_S_N -- 408_S_C:409_S_N:409_S_CA ] 
        [ [ -0.29695226760429244, -0.7603161756971374, -0.5777011889705165, 0.0 ], [ 0.3475891196487567, -0.6495531640853424, 0.6762118683730141, 0.0 ], [ -0.8893824569140494, 0.0, 0.4571639151699629, 15.252466861656293 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.28519, 595, 596, 0, 0, 1,     // 408_S_C:409_S_N:409_S_CA:409_S_C [ 408_S_C:409_S_N:409_S_CA -- 409_S_N:409_S_CA:409_S_C ] 
        [ [ 0.560528359854879, 0.8046881361922669, -0.19566543197462066, -7.692858187848966 ], [ -0.8238990213907825, 0.5657395860410772, -0.033602430527886744, 9.0046586499227 ], [ 0.08365620329485512, 0.18004367319482856, 0.9800948502021647, 21.340211244304307 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 408, ' ') S sidechain

     ],
     [ // 149 : (' ', 409, ' ') S backbone
      [ -80.57654, 596, 597, 0, 0, 1,     // 409_S_N:409_S_CA:409_S_C:409_S_O [ 409_S_N:409_S_CA:409_S_C -- 409_S_CA:409_S_C:409_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 100.39162, 596, 598, 0, 0, 1,     // 409_S_N:409_S_CA:409_S_C:410_M_N [ 409_S_N:409_S_CA:409_S_C -- 409_S_CA:409_S_C:410_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 166.88227, 598, 599, 0, 0, 1,     // 409_S_CA:409_S_C:410_M_N:410_M_CA [ 409_S_CA:409_S_C:410_M_N -- 409_S_C:410_M_N:410_M_CA ] 
        [ [ -0.08030360247880318, -0.983597877520814, -0.1615133021316701, 0.0 ], [ 0.4379008303032195, -0.18037520710050423, 0.880742668140474, 0.0 ], [ -0.8954296143464752, 0.0, 0.4452031061788796, 15.342307820528823 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -68.43095, 599, 600, 0, 0, 1,     // 409_S_C:410_M_N:410_M_CA:410_M_C [ 409_S_C:410_M_N:410_M_CA -- 410_M_N:410_M_CA:410_M_C ] 
        [ [ 0.05733876315712131, 0.9761567820688899, -0.20935664083214894, -2.156010114779348 ], [ -0.9935775088660511, 0.07628569007454852, 0.08357168997562597, 11.756865075302498 ], [ 0.09754998777623583, 0.20322015232441834, 0.9742615509061704, 21.285238758008774 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 409, ' ') S sidechain

     ],
     [ // 150 : (' ', 410, ' ') M backbone
      [ -97.37539, 600, 601, 0, 0, 1,     // 410_M_N:410_M_CA:410_M_C:410_M_O [ 410_M_N:410_M_CA:410_M_C -- 410_M_CA:410_M_C:410_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  84.60739, 600, 602, 0, 0, 1,     // 410_M_N:410_M_CA:410_M_C:411_Y_N [ 410_M_N:410_M_CA:410_M_C -- 410_M_CA:410_M_C:411_Y_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 154.46080, 602, 603, 0, 0, 1,     // 410_M_CA:410_M_C:411_Y_N:411_Y_CA [ 410_M_CA:410_M_C:411_Y_N -- 410_M_C:411_Y_N:411_Y_CA ] 
        [ [ 0.04104521871117962, -0.9955741014192925, 0.08454288026865245, 0.0 ], [ 0.434811984875303, 0.09397983072536507, 0.8956038908053219, 0.0 ], [ -0.8995853643928133, 0.0, 0.43674497383513106, 15.311608692844592 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.42093, 603, 604, 0, 0, 1,     // 410_M_C:411_Y_N:411_Y_CA:411_Y_C [ 410_M_C:411_Y_N:411_Y_CA -- 411_Y_N:411_Y_CA:411_Y_C ] 
        [ [ -0.31877188891176134, 0.8806012958774502, -0.35060781585495643, 1.127674288511346 ], [ -0.9459467866753384, -0.27225697580906183, 0.17624078955198125, 11.946002752006576 ], [ 0.059742444026248794, 0.38783694618004005, 0.9197898366254713, 21.13712677391065 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 410, ' ') M sidechain

     ],
     [ // 151 : (' ', 411, ' ') Y backbone
      [ -58.18656, 604, 605, 0, 0, 1,     // 411_Y_N:411_Y_CA:411_Y_C:411_Y_O [ 411_Y_N:411_Y_CA:411_Y_C -- 411_Y_CA:411_Y_C:411_Y_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 121.06705, 604, 606, 0, 0, 1,     // 411_Y_N:411_Y_CA:411_Y_C:412_P_N [ 411_Y_N:411_Y_CA:411_Y_C -- 411_Y_CA:411_Y_C:412_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.86391, 606, 607, 0, 0, 1,     // 411_Y_CA:411_Y_C:412_P_N:412_P_CA [ 411_Y_CA:411_Y_C:412_P_N -- 411_Y_C:412_P_N:412_P_CA ] 
        [ [ -0.2555933801781897, -0.8565639646327998, -0.44829699809571233, 0.0 ], [ 0.4242533757490585, -0.5160408651381592, 0.7441175301476162, 0.0 ], [ -0.86872383251216, 0.0, 0.4952967825712019, 15.399263081791268 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -75.26483, 607, 608, 0, 0, 1,     // 411_Y_C:412_P_N:412_P_CA:412_P_C [ 411_Y_C:412_P_N:412_P_CA -- 412_P_N:412_P_CA:412_P_C ] 
        [ [ 0.49775656279453534, 0.865495539509413, -0.05617717761048274, -6.025396902099418 ], [ -0.8654862591624547, 0.4998689801769688, 0.03262725636361925, 10.001413081048689 ], [ 0.056319973330515034, 0.03238014431943027, 0.9978875622322902, 22.056366992817598 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 411, ' ') Y sidechain

     ],
     [ // 152 : (' ', 412, ' ') P backbone
      [ -83.91803, 608, 609, 0, 0, 1,     // 412_P_N:412_P_CA:412_P_C:412_P_O [ 412_P_N:412_P_CA:412_P_C -- 412_P_CA:412_P_C:412_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  93.76454, 608, 610, 0, 0, 1,     // 412_P_N:412_P_CA:412_P_C:413_L_N [ 412_P_N:412_P_CA:412_P_C -- 412_P_CA:412_P_C:413_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 171.03156, 610, 611, 0, 0, 1,     // 412_P_CA:412_P_C:413_L_N:413_L_CA [ 412_P_CA:412_P_C:413_L_N -- 412_P_C:413_L_N:413_L_CA ] 
        [ [ -0.02821176846680803, -0.9978422932866134, -0.05928620285097884, 0.0 ], [ 0.42876113376622416, -0.0656563609143257, 0.901028929859075, 0.0 ], [ -0.9029772900197863, 0.0, 0.4296882750419457, 15.298734736405157 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -117.26406, 611, 612, 0, 0, 1,     // 412_P_C:413_L_N:413_L_CA:413_L_C [ 412_P_C:413_L_N:413_L_CA -- 413_L_N:413_L_CA:413_L_C ] 
        [ [ -0.02329324000496028, 0.9900409755184569, -0.1388390858675373, -0.7880349510153377 ], [ -0.988429514026219, -0.001986072918361014, 0.15166789810719733, 11.976518219418656 ], [ 0.14988168924843207, 0.1407654869235532, 0.9786319823709071, 21.010171166659756 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 412, ' ') P sidechain

     ],
     [ // 153 : (' ', 413, ' ') L backbone
      [ -69.70540, 612, 613, 0, 0, 1,     // 413_L_N:413_L_CA:413_L_C:413_L_O [ 413_L_N:413_L_CA:413_L_C -- 413_L_CA:413_L_C:413_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 112.43853, 612, 614, 0, 0, 1,     // 413_L_N:413_L_CA:413_L_C:414_V_N [ 413_L_N:413_L_CA:413_L_C -- 413_L_CA:413_L_C:414_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.32017, 614, 615, 0, 0, 1,     // 413_L_CA:413_L_C:414_V_N:414_V_CA [ 413_L_CA:413_L_C:414_V_N -- 413_L_C:414_V_N:414_V_CA ] 
        [ [ -0.16766771236909056, -0.9242895646556765, -0.34289406366622344, 0.0 ], [ 0.4060171745878541, -0.3816920233232811, 0.8303380355440042, 0.0 ], [ -0.8983527103363199, 0.0, 0.4392748659226795, 15.333988344322188 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  97.54485, 615, 616, 0, 0, 1,     // 413_L_C:414_V_N:414_V_CA:414_V_C [ 413_L_C:414_V_N:414_V_CA -- 414_V_N:414_V_CA:414_V_C ] 
        [ [ 0.3987522155214174, 0.9154394776886764, -0.054472316876048976, -4.603237961942257 ], [ -0.8879365066027111, 0.40025787404570295, 0.22663273043977594, 11.147010028675592 ], [ 0.2292715221280983, -0.04200234461905399, 0.97245584588071, 21.231106296503167 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 413, ' ') L sidechain

     ],
     [ // 154 : (' ', 414, ' ') V backbone
      [ -89.13250, 616, 617, 0, 0, 1,     // 414_V_N:414_V_CA:414_V_C:414_V_O [ 414_V_N:414_V_CA:414_V_C -- 414_V_CA:414_V_C:414_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  98.79092, 616, 618, 0, 0, 1,     // 414_V_N:414_V_CA:414_V_C:415_T_N [ 414_V_N:414_V_CA:414_V_C -- 414_V_CA:414_V_C:415_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 161.47103, 618, 619, 0, 0, 1,     // 414_V_CA:414_V_C:415_T_N:415_T_CA [ 414_V_CA:414_V_C:415_T_N -- 414_V_C:415_T_N:415_T_CA ] 
        [ [ -0.06409804167702438, -0.9882526122207007, -0.13873793818610644, 0.0 ], [ 0.4144825954418225, -0.15282923293323644, 0.8971328795874384, 0.0 ], [ -0.9077971244330867, 0.0, 0.41940956220741915, 15.481783601943006 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  49.54960, 619, 620, 0, 0, 1,     // 414_V_C:415_T_N:415_T_CA:415_T_C [ 414_V_C:415_T_N:415_T_CA -- 415_T_N:415_T_CA:415_T_C ] 
        [ [ -0.04657688983357052, 0.9573939759601873, -0.2850041545812551, -1.8611796236012639 ], [ -0.9788243889358629, 0.013191019452256795, 0.20427631441296584, 12.035103425791556 ], [ 0.199332408197378, 0.28848357284419973, 0.9365061768302861, 21.10819366201514 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 414, ' ') V sidechain

     ],
     [ // 155 : (' ', 415, ' ') T backbone
      [ -87.67675, 620, 621, 0, 0, 1,     // 415_T_N:415_T_CA:415_T_C:415_T_O [ 415_T_N:415_T_CA:415_T_C -- 415_T_CA:415_T_C:415_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 103.16016, 620, 622, 0, 0, 1,     // 415_T_N:415_T_CA:415_T_C:416_A_N [ 415_T_N:415_T_CA:415_T_C -- 415_T_CA:415_T_C:416_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -176.79987, 622, 623, 0, 0, 1,     // 415_T_CA:415_T_C:416_A_N:416_A_CA [ 415_T_CA:415_T_C:416_A_N -- 415_T_C:416_A_N:416_A_CA ] 
        [ [ -0.09402168610277793, -0.9737374292308836, -0.20735318048494145, 0.0 ], [ 0.40212085092252225, -0.2276739311270177, 0.8868277185217638, 0.0 ], [ -0.9107462565367686, 0.0, 0.41296641050364213, 15.41634275879602 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.37492, 623, 624, 0, 0, 1,     // 415_T_C:416_A_N:416_A_CA:416_A_C [ 415_T_C:416_A_N:416_A_CA -- 416_A_N:416_A_CA:416_A_C ] 
        [ [ 0.2545169720221484, 0.9669703618651625, -0.013763365396350685, -2.7769517731276445 ], [ -0.9353860696128887, 0.2497668404086427, 0.2503386230816885, 11.876730319005654 ], [ 0.24550766123856302, -0.05084136806420676, 0.9680605061495541, 20.946944234057842 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 415, ' ') T sidechain

     ],
     [ // 156 : (' ', 416, ' ') A backbone
      [ -165.55215, 624, 625, 0, 0, 1,     // 416_A_N:416_A_CA:416_A_C:416_A_O [ 416_A_N:416_A_CA:416_A_C -- 416_A_CA:416_A_C:416_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  12.12807, 624, 626, 0, 0, 1,     // 416_A_N:416_A_CA:416_A_C:417_T_N [ 416_A_N:416_A_CA:416_A_C -- 416_A_CA:416_A_C:417_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.53417, 626, 627, 0, 0, 1,     // 416_A_CA:416_A_C:417_T_N:417_T_CA [ 416_A_CA:416_A_C:417_T_N -- 416_A_C:417_T_N:417_T_CA ] 
        [ [ 0.4873233447575155, -0.21009764894159313, 0.8475700180891033, 0.0 ], [ 0.10472286065274287, 0.9776804068350838, 0.1821377076486116, 0.0 ], [ -0.8669193042671584, 0.0, 0.49844851277634095, 15.538359727406986 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -106.44582, 627, 628, 0, 0, 1,     // 416_A_C:417_T_N:417_T_CA:417_T_C [ 416_A_C:417_T_N:417_T_CA -- 417_T_N:417_T_CA:417_T_C ] 
        [ [ -0.9843556322952367, 0.1627238344416836, 0.0675643609696736, 11.319151366981199 ], [ -0.15604139968722536, -0.9832101435729791, 0.09459860019500207, 2.4324176628589376 ], [ 0.08182341200595589, 0.08257582745451476, 0.9932200974446227, 22.19505362212771 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 416, ' ') A sidechain

     ],
     [ // 157 : (' ', 417, ' ') T backbone
      [ 148.78793, 628, 629, 0, 0, 1,     // 417_T_N:417_T_CA:417_T_C:417_T_O [ 417_T_N:417_T_CA:417_T_C -- 417_T_CA:417_T_C:417_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -29.67834, 628, 630, 0, 0, 1,     // 417_T_N:417_T_CA:417_T_C:418_Q_N [ 417_T_N:417_T_CA:417_T_C -- 417_T_CA:417_T_C:418_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.77338, 630, 631, 0, 0, 1,     // 417_T_CA:417_T_C:418_Q_N:418_Q_CA [ 417_T_CA:417_T_C:418_Q_N -- 417_T_C:418_Q_N:418_Q_CA ] 
        [ [ 0.4349365589748267, 0.4951302211824558, 0.752114521691311, 0.0 ], [ -0.247865539744695, 0.8688187751607422, -0.42862175648846695, 0.0 ], [ -0.8656748026101997, 0.0, 0.5006067679584365, 15.267850521295852 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -89.29081, 631, 632, 0, 0, 1,     // 417_T_C:418_Q_N:418_Q_CA:418_Q_C [ 417_T_C:418_Q_N:418_Q_CA -- 418_Q_N:418_Q_CA:418_Q_C ] 
        [ [ -0.8430047875031752, -0.5326922496386733, 0.07471208350471623, 10.08176363320925 ], [ 0.5365298321763232, -0.8426269844792553, 0.04599461068694035, -5.7454857104078965 ], [ 0.03845344498966012, 0.07885893863276426, 0.9961438652956415, 21.978263579658815 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 417, ' ') T sidechain

     ],
     [ // 158 : (' ', 418, ' ') Q backbone
      [ 129.73818, 632, 633, 0, 0, 1,     // 418_Q_N:418_Q_CA:418_Q_C:418_Q_O [ 418_Q_N:418_Q_CA:418_Q_C -- 418_Q_CA:418_Q_C:418_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.57310, 632, 634, 0, 0, 1,     // 418_Q_N:418_Q_CA:418_Q_C:419_D_N [ 418_Q_N:418_Q_CA:418_Q_C -- 418_Q_CA:418_Q_C:419_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.00140, 634, 635, 0, 0, 1,     // 418_Q_CA:418_Q_C:419_D_N:419_D_CA [ 418_Q_CA:418_Q_C:419_D_N -- 418_Q_C:419_D_N:419_D_CA ] 
        [ [ 0.34212426532246887, 0.6765303099636382, 0.6521178779776434, 0.0 ], [ -0.3143030790487354, 0.7364147878067793, -0.5990883363912043, 0.0 ], [ -0.8855306666502551, 0.0, 0.464580927742363, 15.371619395499327 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.93730, 635, 636, 0, 0, 1,     // 418_Q_C:419_D_N:419_D_CA:419_D_C [ 418_Q_C:419_D_N:419_D_CA -- 419_D_N:419_D_CA:419_D_C ] 
        [ [ -0.703844980491274, -0.6987404879815372, 0.12792174909905418, 8.742933099839068 ], [ 0.7003476256349622, -0.7127051899521875, -0.039553956599652856, -8.03196695389748 ], [ 0.11880844542669945, 0.06174983943735859, 0.9909952122108115, 21.6002478503872 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 418, ' ') Q sidechain

     ],
     [ // 159 : (' ', 419, ' ') D backbone
      [ 147.25292, 636, 637, 0, 0, 1,     // 419_D_N:419_D_CA:419_D_C:419_D_O [ 419_D_N:419_D_CA:419_D_C -- 419_D_CA:419_D_C:419_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.21536, 636, 638, 0, 0, 1,     // 419_D_N:419_D_CA:419_D_C:420_A_N [ 419_D_N:419_D_CA:419_D_C -- 419_D_CA:419_D_C:420_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -175.27355, 638, 639, 0, 0, 1,     // 419_D_CA:419_D_C:420_A_N:420_A_CA [ 419_D_CA:419_D_C:420_A_N -- 419_D_C:420_A_N:420_A_CA ] 
        [ [ 0.39490985859219896, 0.5477875672603617, 0.7375465983527171, 0.0 ], [ -0.2585730277227279, 0.8366174640509094, -0.4829194634505082, 0.0 ], [ -0.8815816427994578, 0.0, 0.4720315742394879, 15.410517362196837 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -71.45894, 639, 640, 0, 0, 1,     // 419_D_C:420_A_N:420_A_CA:420_A_C [ 419_D_C:420_A_N:420_A_CA -- 420_A_N:420_A_CA:420_A_C ] 
        [ [ -0.8569583306315707, -0.5133847930614673, 0.04537040681309672, 9.871860112552303 ], [ 0.5063515280653264, -0.8550785533197185, -0.11157418016085184, -6.463745340917914 ], [ 0.09607574921415164, -0.07264104835347268, 0.9927198640638999, 21.728531274819616 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 419, ' ') D sidechain

     ],
     [ // 160 : (' ', 420, ' ') A backbone
      [ 154.40579, 640, 641, 0, 0, 1,     // 420_A_N:420_A_CA:420_A_C:420_A_O [ 420_A_N:420_A_CA:420_A_C -- 420_A_CA:420_A_C:420_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -25.61639, 640, 642, 0, 0, 1,     // 420_A_N:420_A_CA:420_A_C:421_D_N [ 420_A_N:420_A_CA:420_A_C -- 420_A_CA:420_A_C:421_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.45854, 642, 643, 0, 0, 1,     // 420_A_CA:420_A_C:421_D_N:421_D_CA [ 420_A_CA:420_A_C:421_D_N -- 420_A_C:421_D_N:421_D_CA ] 
        [ [ 0.4330034311422406, 0.43234363764665035, 0.7909405841183703, 0.0 ], [ -0.2076127607084446, 0.9017089214300045, -0.3792333880403415, 0.0 ], [ -0.8771573235230183, 0.0, 0.48020311305731345, 15.398176386883273 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -85.15140, 643, 644, 0, 0, 1,     // 420_A_C:421_D_N:421_D_CA:421_D_C [ 420_A_C:421_D_N:421_D_CA -- 421_D_N:421_D_CA:421_D_C ] 
        [ [ -0.8841984614136924, -0.4582649631451708, 0.09047819840831503, 10.549489255787785 ], [ 0.4612420899461821, -0.8871625909566115, 0.014080897316223516, -5.0581783675028165 ], [ 0.07381609103533701, 0.054182661070745576, 0.9957988872978095, 21.80307921209451 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 420, ' ') A sidechain

     ],
     [ // 161 : (' ', 421, ' ') D backbone
      [ 126.98748, 644, 645, 0, 0, 1,     // 421_D_N:421_D_CA:421_D_C:421_D_O [ 421_D_N:421_D_CA:421_D_C -- 421_D_CA:421_D_C:421_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.86062, 644, 646, 0, 0, 1,     // 421_D_N:421_D_CA:421_D_C:422_S_N [ 421_D_N:421_D_CA:421_D_C -- 421_D_CA:421_D_C:422_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.55462, 646, 647, 0, 0, 1,     // 421_D_CA:421_D_C:422_S_N:422_S_CA [ 421_D_CA:421_D_C:422_S_N -- 421_D_C:422_S_N:422_S_CA ] 
        [ [ 0.3205666819587662, 0.7296924400631708, 0.6039751197960087, 0.0 ], [ -0.34209339502740993, 0.6837755062260243, -0.6445332932936548, 0.0 ], [ -0.8832944647718379, 0.0, 0.4688186093825973, 15.297475980069805 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.63506, 647, 648, 0, 0, 1,     // 421_D_C:422_S_N:422_S_CA:422_S_C [ 421_D_C:422_S_N:422_S_CA -- 422_S_N:422_S_CA:422_S_C ] 
        [ [ -0.6418949321558565, -0.7568201942803464, 0.12326430790008465, 8.079500252502806 ], [ 0.7612879558539737, -0.6482261062357495, -0.015606519987396325, -8.622055338424161 ], [ 0.09171447183681782, 0.0838218869025115, 0.9922511511365437, 21.56895960370655 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 421, ' ') D sidechain

     ],
     [ // 162 : (' ', 422, ' ') S backbone
      [ 129.64906, 648, 649, 0, 0, 1,     // 422_S_N:422_S_CA:422_S_C:422_S_O [ 422_S_N:422_S_CA:422_S_C -- 422_S_CA:422_S_C:422_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.20027, 648, 650, 0, 0, 1,     // 422_S_N:422_S_CA:422_S_C:423_S_N [ 422_S_N:422_S_CA:422_S_C -- 422_S_CA:422_S_C:423_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.43798, 650, 651, 0, 0, 1,     // 422_S_CA:422_S_C:423_S_N:423_S_CA [ 422_S_CA:422_S_C:423_S_N -- 422_S_C:423_S_N:423_S_CA ] 
        [ [ 0.339268896805378, 0.7095741002607382, 0.6175768874396335, 0.0 ], [ -0.3416490088939819, 0.7046308226576268, -0.6219094455647027, 0.0 ], [ -0.876454545531154, 0.0, 0.48148460994903913, 15.300651640435628 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.48900, 651, 652, 0, 0, 1,     // 422_S_C:423_S_N:423_S_CA:423_S_C [ 422_S_C:423_S_N:423_S_CA -- 423_S_N:423_S_CA:423_S_C ] 
        [ [ -0.6855893017955317, -0.7240304101332905, 0.07581078066960176, 8.27766361072903 ], [ 0.7247249741309635, -0.6886545094799895, -0.02299300858411713, -8.335734855723219 ], [ 0.06885507341067627, 0.03917820535826176, 0.9968570845866103, 21.75420862412505 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 422, ' ') S sidechain

     ],
     [ // 163 : (' ', 423, ' ') S backbone
      [ 142.60919, 652, 653, 0, 0, 1,     // 423_S_N:423_S_CA:423_S_C:423_S_O [ 423_S_N:423_S_CA:423_S_C -- 423_S_CA:423_S_C:423_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -36.81351, 652, 654, 0, 0, 1,     // 423_S_N:423_S_CA:423_S_C:424_R_N [ 423_S_N:423_S_CA:423_S_C -- 423_S_CA:423_S_C:424_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.34349, 654, 655, 0, 0, 1,     // 423_S_CA:423_S_C:424_R_N:424_R_CA [ 423_S_CA:423_S_C:424_R_N -- 423_S_C:424_R_N:424_R_CA ] 
        [ [ 0.3765179820289938, 0.5992123310640113, 0.7065258604677171, 0.0 ], [ -0.2818098858085228, 0.8005901462676354, -0.5288086666836632, 0.0 ], [ -0.8825063158240861, 0.0, 0.47030054489719475, 15.293371296586947 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.87087, 655, 656, 0, 0, 1,     // 423_S_C:424_R_N:424_R_CA:424_R_C [ 423_S_C:424_R_N:424_R_CA -- 424_R_N:424_R_CA:424_R_C ] 
        [ [ -0.7717515824583457, -0.6278010052764492, 0.1013182744963558, 9.493107605544735 ], [ 0.6317983445746906, -0.7750697275507737, 0.009887832272522356, -7.105242506267795 ], [ 0.07232113636911056, 0.07164366830536352, 0.994804924608856, 21.61247995129328 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 423, ' ') S sidechain

     ],
     [ // 164 : (' ', 424, ' ') R backbone
      [ 128.35316, 656, 657, 0, 0, 1,     // 424_R_N:424_R_CA:424_R_C:424_R_O [ 424_R_N:424_R_CA:424_R_C -- 424_R_CA:424_R_C:424_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.59738, 656, 658, 0, 0, 1,     // 424_R_N:424_R_CA:424_R_C:425_K_N [ 424_R_N:424_R_CA:424_R_C -- 424_R_CA:424_R_C:425_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.39201, 658, 659, 0, 0, 1,     // 424_R_CA:424_R_C:425_K_N:425_K_CA [ 424_R_CA:424_R_C:425_K_N -- 424_R_C:425_K_N:425_K_CA ] 
        [ [ 0.3151271290455028, 0.7500808184523089, 0.5814410187881944, 0.0 ], [ -0.3574086055123418, 0.6613461769677923, -0.6594545647092064, 0.0 ], [ -0.8791780145370233, 0.0, 0.4764934613976752, 15.199180382819295 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -55.46948, 659, 660, 0, 0, 1,     // 424_R_C:425_K_N:425_K_CA:425_K_C [ 424_R_C:425_K_N:425_K_CA -- 425_K_N:425_K_CA:425_K_C ] 
        [ [ -0.6332112407068555, -0.7684250053483296, 0.09255558220822985, 7.781046891727404 ], [ 0.7700027197630442, -0.6375437180961511, -0.02517576361602947, -8.825051424236054 ], [ 0.07835391630288456, 0.055326473514540474, 0.9953891927926711, 21.575782388752785 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 424, ' ') R sidechain

     ],
     [ // 165 : (' ', 425, ' ') K backbone
      [ 122.36934, 660, 661, 0, 0, 1,     // 425_K_N:425_K_CA:425_K_C:425_K_O [ 425_K_N:425_K_CA:425_K_C -- 425_K_CA:425_K_C:425_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -51.24703, 660, 662, 0, 0, 1,     // 425_K_N:425_K_CA:425_K_C:426_L_N [ 425_K_N:425_K_CA:425_K_C -- 425_K_CA:425_K_C:426_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.90334, 662, 663, 0, 0, 1,     // 425_K_CA:425_K_C:426_L_N:426_L_CA [ 425_K_CA:425_K_C:426_L_N -- 425_K_C:426_L_N:426_L_CA ] 
        [ [ 0.29726754421591617, 0.7798520764347117, 0.5508745283967186, 0.0 ], [ -0.370348403330902, 0.6259638479021522, -0.6863026455360514, 0.0 ], [ -0.8800420826904188, 0.0, 0.474895707175702, 15.242002569124724 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.81126, 663, 664, 0, 0, 1,     // 425_K_C:426_L_N:426_L_CA:426_L_C [ 425_K_C:426_L_N:426_L_CA -- 426_L_N:426_L_CA:426_L_C ] 
        [ [ -0.6076852342470156, -0.790205633641404, 0.07933292279714156, 7.347036669489281 ], [ 0.7892755354669184, -0.6119954502793583, -0.05005694707805087, -9.153247142858655 ], [ 0.08810666939317963, 0.032196667510053796, 0.9955905731824142, 21.575706516115407 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 425, ' ') K sidechain

     ],
     [ // 166 : (' ', 426, ' ') L backbone
      [ 133.92667, 664, 665, 0, 0, 1,     // 426_L_N:426_L_CA:426_L_C:426_L_O [ 426_L_N:426_L_CA:426_L_C -- 426_L_CA:426_L_C:426_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.42544, 664, 666, 0, 0, 1,     // 426_L_N:426_L_CA:426_L_C:427_A_N [ 426_L_N:426_L_CA:426_L_C -- 426_L_CA:426_L_C:427_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.78231, 666, 667, 0, 0, 1,     // 426_L_CA:426_L_C:427_A_N:427_A_CA [ 426_L_CA:426_L_C:427_A_N -- 426_L_C:427_A_N:427_A_CA ] 
        [ [ 0.3500451151808566, 0.6746302541722279, 0.6498787867699145, 0.0 ], [ -0.31992028975980563, 0.7381558237631911, -0.5939503245596275, 0.0 ], [ -0.8804086696177081, 0.0, 0.4742157467461168, 15.3947879637408 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.28441, 667, 668, 0, 0, 1,     // 426_L_C:427_A_N:427_A_CA:427_A_C [ 426_L_C:427_A_N:427_A_CA -- 427_A_N:427_A_CA:427_A_C ] 
        [ [ -0.7023783914650322, -0.7036679441890568, 0.10731271836898201, 8.690767838592924 ], [ 0.70818454448722, -0.7060036262975701, 0.005790561640264891, -7.942841778327152 ], [ 0.07168853571124596, 0.08006437394639122, 0.9942084539331527, 21.736430504282133 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 426, ' ') L sidechain

     ],
     [ // 167 : (' ', 427, ' ') A backbone
      [ 139.61289, 668, 669, 0, 0, 1,     // 427_A_N:427_A_CA:427_A_C:427_A_O [ 427_A_N:427_A_CA:427_A_C -- 427_A_CA:427_A_C:427_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.97530, 668, 670, 0, 0, 1,     // 427_A_N:427_A_CA:427_A_C:428_H_N [ 427_A_N:427_A_CA:427_A_C -- 427_A_CA:427_A_C:428_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.01307, 670, 671, 0, 0, 1,     // 427_A_CA:427_A_C:428_H_N:428_H_CA [ 427_A_CA:427_A_C:428_H_N -- 427_A_C:428_H_N:428_H_CA ] 
        [ [ 0.3834320826054336, 0.6153216992733637, 0.688737282591991, 0.0 ], [ -0.2993038622831154, 0.7882760978257175, -0.5376225363757268, 0.0 ], [ -0.8737259501990712, 0.0, 0.4864185069965268, 15.42165803975486 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.15744, 671, 672, 0, 0, 1,     // 427_A_C:428_H_N:428_H_CA:428_H_C [ 427_A_C:428_H_N:428_H_CA -- 428_H_N:428_H_CA:428_H_C ] 
        [ [ -0.7446362779220366, -0.6573940691244397, 0.1155415573827316, 9.217657304990544 ], [ 0.6647525985272678, -0.7460141114778163, 0.03958444425793364, -7.195225850270951 ], [ 0.06017305338488909, 0.1062825637438983, 0.9925135869550419, 21.931599083494238 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 427, ' ') A sidechain

     ],
     [ // 168 : (' ', 428, ' ') H backbone
      [ 133.13577, 672, 673, 0, 0, 1,     // 428_H_N:428_H_CA:428_H_C:428_H_O [ 428_H_N:428_H_CA:428_H_C -- 428_H_CA:428_H_C:428_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.40160, 672, 674, 0, 0, 1,     // 428_H_N:428_H_CA:428_H_C:429_L_N [ 428_H_N:428_H_CA:428_H_C -- 428_H_CA:428_H_C:429_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.80409, 674, 675, 0, 0, 1,     // 428_H_CA:428_H_C:429_L_N:429_L_CA [ 428_H_CA:428_H_C:429_L_N -- 428_H_C:429_L_N:429_L_CA ] 
        [ [ 0.35882292321402365, 0.6743230654437636, 0.6453948506043944, 0.0 ], [ -0.32766878004275485, 0.7384364586140947, -0.5893596246560898, 0.0 ], [ -0.8740018766349622, 0.0, 0.48592254489636966, 15.394623031389012 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.11210, 675, 676, 0, 0, 1,     // 428_H_C:429_L_N:429_L_CA:429_L_C [ 428_H_C:429_L_N:429_L_CA -- 429_L_N:429_L_CA:429_L_C ] 
        [ [ -0.7031337793303474, -0.7040477432018324, 0.09959750828725097, 8.664933839221765 ], [ 0.7083774158442607, -0.7057279203851936, 0.01222861850729652, -7.912616827313156 ], [ 0.06167921113657045, 0.07915098029207464, 0.9949526607996855, 21.91851565087661 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 428, ' ') H sidechain

     ],
     [ // 169 : (' ', 429, ' ') L backbone
      [ 131.08800, 676, 677, 0, 0, 1,     // 429_L_N:429_L_CA:429_L_C:429_L_O [ 429_L_N:429_L_CA:429_L_C -- 429_L_CA:429_L_C:429_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -47.53245, 676, 678, 0, 0, 1,     // 429_L_N:429_L_CA:429_L_C:430_L_N [ 429_L_N:429_L_CA:429_L_C -- 429_L_CA:429_L_C:430_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.32664, 678, 679, 0, 0, 1,     // 429_L_CA:429_L_C:430_L_N:430_L_CA [ 429_L_CA:429_L_C:430_L_N -- 429_L_C:430_L_N:430_L_CA ] 
        [ [ 0.309982889856392, 0.7376598067808252, 0.5998071502210987, 0.0 ], [ -0.33867180961413607, 0.6751725775385695, -0.6553193083624278, 0.0 ], [ -0.8883760540272156, 0.0, 0.4591165283792707, 15.372894672731622 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.67326, 679, 680, 0, 0, 1,     // 429_L_C:430_L_N:430_L_CA:430_L_C [ 429_L_C:430_L_N:430_L_CA -- 430_L_N:430_L_CA:430_L_C ] 
        [ [ -0.676370240018912, -0.7339659339768941, 0.06178435221142107, 8.025378867116379 ], [ 0.7294380697836049, -0.6791060368696548, -0.08206761259814947, -8.768127766743907 ], [ 0.10219295850071174, -0.010440232210860052, 0.9947098073228469, 21.515842590309195 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 429, ' ') L sidechain

     ],
     [ // 170 : (' ', 430, ' ') L backbone
      [ 138.71848, 680, 681, 0, 0, 1,     // 430_L_N:430_L_CA:430_L_C:430_L_O [ 430_L_N:430_L_CA:430_L_C -- 430_L_CA:430_L_C:430_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.57450, 680, 682, 0, 0, 1,     // 430_L_N:430_L_CA:430_L_C:431_N_N [ 430_L_N:430_L_CA:430_L_C -- 430_L_CA:430_L_C:431_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.23389, 682, 683, 0, 0, 1,     // 430_L_CA:430_L_C:431_N_N:431_N_CA [ 430_L_CA:430_L_C:431_N_N -- 430_L_C:431_N_N:431_N_CA ] 
        [ [ 0.3700820014614894, 0.6235316543888149, 0.6886563643570034, 0.0 ], [ -0.2951629600340288, 0.7817981043563278, -0.5492454378953027, 0.0 ], [ -0.8808621567635929, 0.0, 0.47337285598330575, 15.37429927240189 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.24732, 683, 684, 0, 0, 1,     // 430_L_C:431_N_N:431_N_CA:431_N_C [ 430_L_C:431_N_N:431_N_CA -- 431_N_N:431_N_CA:431_N_C ] 
        [ [ -0.7374927159935644, -0.6627908020017623, 0.12966436148108199, 9.207367530230997 ], [ 0.6704777314253209, -0.7415780602080126, 0.02283839488397715, -7.343436977782845 ], [ 0.08101916760366716, 0.10378021680446645, 0.9912948910797568, 21.703316311172657 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 430, ' ') L sidechain

     ],
     [ // 171 : (' ', 431, ' ') N backbone
      [ 126.67006, 684, 685, 0, 0, 1,     // 431_N_N:431_N_CA:431_N_C:431_N_O [ 431_N_N:431_N_CA:431_N_C -- 431_N_CA:431_N_C:431_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -47.80956, 684, 686, 0, 0, 1,     // 431_N_N:431_N_CA:431_N_C:432_A_N [ 431_N_N:431_N_CA:431_N_C -- 431_N_CA:431_N_C:432_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.28776, 686, 687, 0, 0, 1,     // 431_N_CA:431_N_C:432_A_N:432_A_CA [ 431_N_CA:431_N_C:432_A_N -- 431_N_C:432_A_N:432_A_CA ] 
        [ [ 0.3284025191874322, 0.7409167232833462, 0.585827785744597, 0.0 ], [ -0.3622990441277041, 0.6715969097301366, -0.6462948193085587, 0.0 ], [ -0.8722907703372197, 0.0, 0.48898753765765834, 15.346965803918998 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.17770, 687, 688, 0, 0, 1,     // 431_N_C:432_A_N:432_A_CA:432_A_C [ 431_N_C:432_A_N:432_A_CA -- 432_A_N:432_A_CA:432_A_C ] 
        [ [ -0.66565665704574, -0.7449417004675791, 0.04430663421161577, 7.831282122590418 ], [ 0.7442782498553365, -0.6670414463830626, -0.03325049773265395, -8.639598850609143 ], [ 0.05432404369123808, 0.010843049002176577, 0.9984644843785726, 21.883698022860607 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 431, ' ') N sidechain

     ],
     [ // 172 : (' ', 432, ' ') A backbone
      [ 137.04562, 688, 689, 0, 0, 1,     // 432_A_N:432_A_CA:432_A_C:432_A_O [ 432_A_N:432_A_CA:432_A_C -- 432_A_CA:432_A_C:432_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.36343, 688, 690, 0, 0, 1,     // 432_A_N:432_A_CA:432_A_C:433_V_N [ 432_A_N:432_A_CA:432_A_C -- 432_A_CA:432_A_C:433_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.20505, 690, 691, 0, 0, 1,     // 432_A_CA:432_A_C:433_V_N:433_V_CA [ 432_A_CA:432_A_C:433_V_N -- 432_A_C:433_V_N:433_V_CA ] 
        [ [ 0.3694399639667784, 0.6342371932547964, 0.6791592565197824, 0.0 ], [ -0.3030667305929988, 0.7731385274918577, -0.5571421489309463, 0.0 ], [ -0.8784444603000785, 0.0, 0.47784446231813116, 15.416608807595585 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.00408, 691, 692, 0, 0, 1,     // 432_A_C:433_V_N:433_V_CA:433_V_C [ 432_A_C:433_V_N:433_V_CA -- 433_V_N:433_V_CA:433_V_C ] 
        [ [ -0.7297362187446398, -0.6734931050608762, 0.11786470416428518, 9.098936587649304 ], [ 0.680690897875035, -0.731850167352984, 0.03248436693992817, -7.464230274069281 ], [ 0.06438130631133296, 0.10393445040441922, 0.9924982002083317, 21.818459948772528 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 432, ' ') A sidechain

     ],
     [ // 173 : (' ', 433, ' ') V backbone
      [ 144.27228, 692, 693, 0, 0, 1,     // 433_V_N:433_V_CA:433_V_C:433_V_O [ 433_V_N:433_V_CA:433_V_C -- 433_V_CA:433_V_C:433_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -32.91453, 692, 694, 0, 0, 1,     // 433_V_N:433_V_CA:433_V_C:434_T_N [ 433_V_N:433_V_CA:433_V_C -- 433_V_CA:433_V_C:434_T_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.09809, 694, 695, 0, 0, 1,     // 433_V_CA:433_V_C:434_T_N:434_T_CA [ 433_V_CA:433_V_C:434_T_N -- 433_V_C:434_T_N:434_T_CA ] 
        [ [ 0.4015940484196004, 0.5433873152403407, 0.7371922720090399, 0.0 ], [ -0.25994730231048935, 0.8394821175200188, -0.47717625086083704, 0.0 ], [ -0.878151251377264, 0.0, 0.4783830888571891, 15.385351613827325 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -68.19118, 695, 696, 0, 0, 1,     // 433_V_C:434_T_N:434_T_CA:434_T_C [ 433_V_C:434_T_N:434_T_CA -- 434_T_N:434_T_CA:434_T_C ] 
        [ [ -0.8049789882196836, -0.5818012260674026, 0.11625902920324677, 9.841430514802804 ], [ 0.5881365534401579, -0.8083030802192169, 0.02723095693614563, -6.37024707728278 ], [ 0.07812952727584209, 0.09029653290460853, 0.9928455635762595, 21.771709837054914 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 433, ' ') V sidechain

     ],
     [ // 174 : (' ', 434, ' ') T backbone
      [ 133.05048, 696, 697, 0, 0, 1,     // 434_T_N:434_T_CA:434_T_C:434_T_O [ 434_T_N:434_T_CA:434_T_C -- 434_T_CA:434_T_C:434_T_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.69139, 696, 698, 0, 0, 1,     // 434_T_N:434_T_CA:434_T_C:435_D_N [ 434_T_N:434_T_CA:434_T_C -- 434_T_CA:434_T_C:435_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.59162, 698, 699, 0, 0, 1,     // 434_T_CA:434_T_C:435_D_N:435_D_CA [ 434_T_CA:434_T_C:435_D_N -- 434_T_C:435_D_N:435_D_CA ] 
        [ [ 0.3509351587789224, 0.6780492022369041, 0.6458279907829254, 0.0 ], [ -0.32373599660054386, 0.7350165163762634, -0.5957732162149914, 0.0 ], [ -0.8786577939322369, 0.0, 0.4774520721100023, 15.32308070804763 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.48935, 699, 700, 0, 0, 1,     // 434_T_C:435_D_N:435_D_CA:435_D_C [ 434_T_C:435_D_N:435_D_CA -- 435_D_N:435_D_CA:435_D_C ] 
        [ [ -0.7008632645922004, -0.7030178119929333, 0.12065090288030086, 8.611837332751705 ], [ 0.7058069883659632, -0.7079580340313991, -0.02513796380308646, -7.944378531865674 ], [ 0.10308821231805448, 0.06753797502934773, 0.9923766635758838, 21.689697884568563 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 434, ' ') T sidechain

     ],
     [ // 175 : (' ', 435, ' ') D backbone
      [ 137.70859, 700, 701, 0, 0, 1,     // 435_D_N:435_D_CA:435_D_C:435_D_O [ 435_D_N:435_D_CA:435_D_C -- 435_D_CA:435_D_C:435_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.18088, 700, 702, 0, 0, 1,     // 435_D_N:435_D_CA:435_D_C:436_A_N [ 435_D_N:435_D_CA:435_D_C -- 435_D_CA:435_D_C:436_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 168.80958, 702, 703, 0, 0, 1,     // 435_D_CA:435_D_C:436_A_N:436_A_CA [ 435_D_CA:435_D_C:436_A_N -- 435_D_C:436_A_N:436_A_CA ] 
        [ [ 0.3742995715052142, 0.6043332693316933, 0.7033357166743915, 0.0 ], [ -0.2839120145056906, 0.7967316358591938, -0.5334910200185693, 0.0 ], [ -0.8827761883911083, 0.0, 0.46979378583551573, 15.387607326150572 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.74856, 703, 704, 0, 0, 1,     // 435_D_C:436_A_N:436_A_CA:436_A_C [ 435_D_C:436_A_N:436_A_CA -- 436_A_N:436_A_CA:436_A_C ] 
        [ [ -0.7230509555854445, -0.6654839021795931, 0.1852794958078124, 9.416780155763409 ], [ 0.6844193394644161, -0.7264849604296217, 0.06156110815005806, -7.142773403209081 ], [ 0.09363484070615827, 0.1713206882118814, 0.9807556976110512, 21.677554921410987 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 435, ' ') D sidechain

     ],
     [ // 176 : (' ', 436, ' ') A backbone
      [ 130.27115, 704, 705, 0, 0, 1,     // 436_A_N:436_A_CA:436_A_C:436_A_O [ 436_A_N:436_A_CA:436_A_C -- 436_A_CA:436_A_C:436_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.55479, 704, 706, 0, 0, 1,     // 436_A_N:436_A_CA:436_A_C:437_L_N [ 436_A_N:436_A_CA:436_A_C -- 436_A_CA:436_A_C:437_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.66692, 706, 707, 0, 0, 1,     // 436_A_CA:436_A_C:437_L_N:437_L_CA [ 436_A_CA:436_A_C:437_L_N -- 436_A_C:437_L_N:437_L_CA ] 
        [ [ 0.336125013039007, 0.7139204306788794, 0.6142781082448002, 0.0 ], [ -0.34269825507323926, 0.700226833718391, -0.626290896716912, 0.0 ], [ -0.877255881473179, 0.0, 0.4800230394686442, 15.363968768678758 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.37615, 707, 708, 0, 0, 1,     // 436_A_C:437_L_N:437_L_CA:437_L_C [ 436_A_C:437_L_N:437_L_CA -- 437_L_N:437_L_CA:437_L_C ] 
        [ [ -0.6691516109921661, -0.737275559857518, 0.09306379716824544, 8.208618847552016 ], [ 0.7402385891174206, -0.6723329765569376, -0.0038986940812392617, -8.369146140548978 ], [ 0.06544427162127642, 0.06628059648850493, 0.9956525146059184, 21.778532654571354 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 436, ' ') A sidechain

     ],
     [ // 177 : (' ', 437, ' ') L backbone
      [ 131.48146, 708, 709, 0, 0, 1,     // 437_L_N:437_L_CA:437_L_C:437_L_O [ 437_L_N:437_L_CA:437_L_C -- 437_L_CA:437_L_C:437_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.91872, 708, 710, 0, 0, 1,     // 437_L_N:437_L_CA:437_L_C:438_V_N [ 437_L_N:437_L_CA:437_L_C -- 437_L_CA:437_L_C:438_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.34300, 710, 711, 0, 0, 1,     // 437_L_CA:437_L_C:438_V_N:438_V_CA [ 437_L_CA:437_L_C:438_V_N -- 437_L_C:438_V_N:438_V_CA ] 
        [ [ 0.3284345879809659, 0.7183536548068215, 0.6132689035353539, 0.0 ], [ -0.33913987624216907, 0.6956781056104051, -0.6332583341076249, 0.0 ], [ -0.8815411877843369, 0.0, 0.4721071215728276, 15.379415783143829 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.36683, 711, 712, 0, 0, 1,     // 437_L_C:438_V_N:438_V_CA:438_V_C [ 437_L_C:438_V_N:438_V_CA -- 438_V_N:438_V_CA:438_V_C ] 
        [ [ -0.6668546774520774, -0.7378395921410829, 0.10439145285590412, 8.188764010762101 ], [ 0.7393117274866733, -0.6726300491797339, -0.031416341943958384, -8.455675847840768 ], [ 0.09339704899490571, 0.05622769077197956, 0.9940399579642127, 21.683296159392505 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 437, ' ') L sidechain

     ],
     [ // 178 : (' ', 438, ' ') V backbone
      [ 135.47278, 712, 713, 0, 0, 1,     // 438_V_N:438_V_CA:438_V_C:438_V_O [ 438_V_N:438_V_CA:438_V_C -- 438_V_CA:438_V_C:438_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.55229, 712, 714, 0, 0, 1,     // 438_V_N:438_V_CA:438_V_C:439_W_N [ 438_V_N:438_V_CA:438_V_C -- 438_V_CA:438_V_C:439_W_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.02803, 714, 715, 0, 0, 1,     // 438_V_CA:438_V_C:439_W_N:439_W_CA [ 438_V_CA:438_V_C:439_W_N -- 438_V_C:439_W_N:439_W_CA ] 
        [ [ 0.3616804037072748, 0.6501417564764788, 0.6682087862785266, 0.0 ], [ -0.3094755591253902, 0.7598129351923926, -0.5717597238511543, 0.0 ], [ -0.8794385503707294, 0.0, 0.4760124327387155, 15.319401948519312 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.81039, 715, 716, 0, 0, 1,     // 438_V_C:439_W_N:439_W_CA:439_W_C [ 438_V_C:439_W_N:439_W_CA -- 439_W_N:439_W_CA:439_W_C ] 
        [ [ -0.724525131969832, -0.6790416505475457, 0.11817685884200076, 8.921202066670462 ], [ 0.6832773303514945, -0.7301322006067302, -0.0062497571879291305, -7.6335183475628705 ], [ 0.09052857544351373, 0.07621946246752392, 0.9929728951836146, 21.674605869677226 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 438, ' ') V sidechain

     ],
     [ // 179 : (' ', 439, ' ') W backbone
      [ 132.05867, 716, 717, 0, 0, 1,     // 439_W_N:439_W_CA:439_W_C:439_W_O [ 439_W_N:439_W_CA:439_W_C -- 439_W_CA:439_W_C:439_W_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.95716, 716, 718, 0, 0, 1,     // 439_W_N:439_W_CA:439_W_C:440_V_N [ 439_W_N:439_W_CA:439_W_C -- 439_W_CA:439_W_C:440_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.13686, 718, 719, 0, 0, 1,     // 439_W_CA:439_W_C:440_V_N:440_V_CA [ 439_W_CA:439_W_C:440_V_N -- 439_W_C:440_V_N:440_V_CA ] 
        [ [ 0.3504671598649605, 0.6941202681224856, 0.6287844012360333, 0.0 ], [ -0.3379360982098707, 0.7198590510520575, -0.6063020205682312, 0.0 ], [ -0.8734826634701325, 0.0, 0.4868552522230023, 15.399954329968104 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.48010, 719, 720, 0, 0, 1,     // 439_W_C:440_V_N:440_V_CA:440_V_C [ 439_W_C:440_V_N:440_V_CA -- 440_V_N:440_V_CA:440_V_C ] 
        [ [ -0.6847566682342422, -0.7213326255794085, 0.10386312427159493, 8.435377365577082 ], [ 0.7251212271865384, -0.6886187648679141, -0.001844601658060648, -8.133767839900768 ], [ 0.07285266770838539, 0.07405025284565389, 0.9945898897843605, 21.931299202457428 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 439, ' ') W sidechain

     ],
     [ // 180 : (' ', 440, ' ') V backbone
      [ 131.77394, 720, 721, 0, 0, 1,     // 440_V_N:440_V_CA:440_V_C:440_V_O [ 440_V_N:440_V_CA:440_V_C -- 440_V_CA:440_V_C:440_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.25252, 720, 722, 0, 0, 1,     // 440_V_N:440_V_CA:440_V_C:441_I_N [ 440_V_N:440_V_CA:440_V_C -- 440_V_CA:440_V_C:441_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.12121, 722, 723, 0, 0, 1,     // 440_V_CA:440_V_C:441_I_N:441_I_CA [ 440_V_CA:440_V_C:441_I_N -- 440_V_C:441_I_N:441_I_CA ] 
        [ [ 0.3240480654506063, 0.7102163804469401, 0.6249684345809527, 0.0 ], [ -0.326917124656332, 0.7039834464991681, -0.6305017848203042, 0.0 ], [ -0.887760128009901, 0.0, 0.4603063709268471, 15.332388829752615 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.30459, 723, 724, 0, 0, 1,     // 440_V_C:441_I_N:441_I_CA:441_I_C [ 440_V_C:441_I_N:441_I_CA -- 441_I_N:441_I_CA:441_I_C ] 
        [ [ -0.6868046915521229, -0.7204585868859291, 0.0961183658013566, 8.342997678818072 ], [ 0.7186933294704267, -0.6928869291963711, -0.058203105789431546, -8.416864974586296 ], [ 0.10853208666890156, 0.029105462221947268, 0.9936667742418184, 21.477235147891605 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 440, ' ') V sidechain

     ],
     [ // 181 : (' ', 441, ' ') I backbone
      [ 131.86421, 724, 725, 0, 0, 1,     // 441_I_N:441_I_CA:441_I_C:441_I_O [ 441_I_N:441_I_CA:441_I_C -- 441_I_CA:441_I_C:441_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.78682, 724, 726, 0, 0, 1,     // 441_I_N:441_I_CA:441_I_C:442_A_N [ 441_I_N:441_I_CA:441_I_C -- 441_I_CA:441_I_C:442_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.65094, 726, 727, 0, 0, 1,     // 441_I_CA:441_I_C:442_A_N:442_A_CA [ 441_I_CA:441_I_C:442_A_N -- 441_I_C:442_A_N:442_A_CA ] 
        [ [ 0.3395307608796542, 0.7044710007991632, 0.6232491247081768, 0.0 ], [ -0.33701356908853636, 0.7097327729737619, -0.6186285195633814, 0.0 ], [ -0.8781461818323242, 0.0, 0.4783923947277074, 15.383251554927647 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.23552, 727, 728, 0, 0, 1,     // 441_I_C:442_A_N:442_A_CA:442_A_C [ 441_I_C:442_A_N:442_A_CA -- 442_A_N:442_A_CA:442_A_C ] 
        [ [ -0.6911499703199443, -0.7177954600872526, 0.08415103091982493, 8.326936537074358 ], [ 0.7178287564640692, -0.6953230852479217, -0.0353225638163364, -8.265202818921756 ], [ 0.08386653039197446, 0.03599284094706213, 0.9958267522418603, 21.774825891404824 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 441, ' ') I sidechain

     ],
     [ // 182 : (' ', 442, ' ') A backbone
      [ 139.90660, 728, 729, 0, 0, 1,     // 442_A_N:442_A_CA:442_A_C:442_A_O [ 442_A_N:442_A_CA:442_A_C -- 442_A_CA:442_A_C:442_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -36.21824, 728, 730, 0, 0, 1,     // 442_A_N:442_A_CA:442_A_C:443_K_N [ 442_A_N:442_A_CA:442_A_C -- 442_A_CA:442_A_C:443_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.88648, 730, 731, 0, 0, 1,     // 442_A_CA:442_A_C:443_K_N:443_K_CA [ 442_A_CA:442_A_C:443_K_N -- 442_A_C:443_K_N:443_K_CA ] 
        [ [ 0.3816042065384119, 0.5908625843330132, 0.7108161759467083, 0.0 ], [ -0.27947869740089, 0.8067722147144836, -0.5205864493662785, 0.0 ], [ -0.8810617953647126, 0.0, 0.4730011762656719, 15.421593555387377 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -67.56892, 731, 732, 0, 0, 1,     // 442_A_C:443_K_N:443_K_CA:443_K_C [ 442_A_C:443_K_N:443_K_CA -- 443_K_N:443_K_CA:443_K_C ] 
        [ [ -0.7723266935166562, -0.6225231317201512, 0.12639789932971324, 9.537016271784358 ], [ 0.627455777095919, -0.778651543536777, -0.00101071101265382, -6.984705197323574 ], [ 0.09904911039774866, 0.07852849305271316, 0.9919791073950522, 21.76784743545813 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 442, ' ') A sidechain

     ],
     [ // 183 : (' ', 443, ' ') K backbone
      [ 163.41534, 732, 733, 0, 0, 1,     // 443_K_N:443_K_CA:443_K_C:443_K_O [ 443_K_N:443_K_CA:443_K_C -- 443_K_CA:443_K_C:443_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -12.91218, 732, 734, 0, 0, 1,     // 443_K_N:443_K_CA:443_K_C:444_S_N [ 443_K_N:443_K_CA:443_K_C -- 443_K_CA:443_K_C:444_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.55522, 734, 735, 0, 0, 1,     // 443_K_CA:443_K_C:444_S_N:444_S_CA [ 443_K_CA:443_K_C:444_S_N -- 443_K_C:444_S_N:444_S_CA ] 
        [ [ 0.47690279877463587, 0.22345740665959582, 0.850076765904049, 0.0 ], [ -0.10933206661502391, 0.9747136950966515, -0.1948838418153746, 0.0 ], [ -0.8721297035020692, 0.0, 0.48927474926608755, 15.330672055801783 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -77.18107, 735, 736, 0, 0, 1,     // 443_K_C:444_S_N:444_S_CA:444_S_C [ 443_K_C:444_S_N:444_S_CA -- 444_S_N:444_S_CA:444_S_C ] 
        [ [ -0.9593042881863124, -0.26770076178978536, 0.08984199911244643, 11.348896081203787 ], [ 0.27355818102205, -0.9599417732428855, 0.06064415535904271, -2.601784400396243 ], [ 0.07000860135197376, 0.08275321214595209, 0.9941079929345038, 21.862703640822918 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 443, ' ') K sidechain

     ],
     [ // 184 : (' ', 444, ' ') S backbone
      [ 179.22356, 736, 737, 0, 0, 1,     // 444_S_N:444_S_CA:444_S_C:444_S_O [ 444_S_N:444_S_CA:444_S_C -- 444_S_CA:444_S_C:444_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -1.79704, 736, 738, 0, 0, 1,     // 444_S_N:444_S_CA:444_S_C:445_G_N [ 444_S_N:444_S_CA:444_S_C -- 444_S_CA:444_S_C:445_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.61151, 738, 739, 0, 0, 1,     // 444_S_CA:444_S_C:445_G_N:445_G_CA [ 444_S_CA:444_S_C:445_G_N -- 444_S_C:445_G_N:445_G_CA ] 
        [ [ 0.4886023255402184, 0.031359112484115916, 0.8719428728654757, 0.0 ], [ -0.015329674695190868, 0.9995081820896758, -0.02735680919864485, 0.0 ], [ -0.8723719210006877, 0.0, 0.48884274715860304, 15.393901954387616 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  74.78407, 739, 740, 0, 0, 1,     // 444_S_C:445_G_N:445_G_CA:445_G_C [ 444_S_C:445_G_N:445_G_CA -- 445_G_N:445_G_CA:445_G_C ] 
        [ [ -0.992979709912944, -0.0855308564881341, 0.08170537491261805, 11.691919889481532 ], [ 0.09282670897574878, -0.991595795659774, 0.09011648090438919, -0.366828643866597 ], [ 0.07331096645070864, 0.09706827812558379, 0.9925740534488081, 21.94881560320215 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 444, ' ') S sidechain

     ],
     [ // 185 : (' ', 445, ' ') G backbone
      [ -172.39288, 740, 741, 0, 0, 1,     // 445_G_N:445_G_CA:445_G_C:445_G_O [ 445_G_N:445_G_CA:445_G_C -- 445_G_CA:445_G_C:445_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [   8.08127, 740, 742, 0, 0, 1,     // 445_G_N:445_G_CA:445_G_C:446_I_N [ 445_G_N:445_G_CA:445_G_C -- 445_G_CA:445_G_C:446_I_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.91096, 742, 743, 0, 0, 1,     // 445_G_CA:445_G_C:446_I_N:446_I_CA [ 445_G_CA:445_G_C:446_I_N -- 445_G_C:446_I_N:446_I_CA ] 
        [ [ 0.4695350644971696, -0.14057752096265108, 0.8716506087920974, 0.0 ], [ 0.06666811141090358, 0.9900696746189105, 0.12376349349023706, 0.0 ], [ -0.880393199728702, 0.0, 0.47424446635828854, 15.354906185116361 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -80.87721, 743, 744, 0, 0, 1,     // 445_G_C:446_I_N:446_I_CA:446_I_C [ 445_G_C:446_I_N:446_I_CA -- 446_I_N:446_I_CA:446_I_C ] 
        [ [ -0.9866241778225189, 0.13162804874182552, 0.09616022317184096, 11.643418782708245 ], [ -0.129351209312857, -0.9911579418804535, 0.029566854695843463, 1.6532199601337334 ], [ 0.09920179628081394, 0.016732932550012593, 0.9949266367843091, 21.689814772908775 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 445, ' ') G sidechain

     ],
     [ // 186 : (' ', 446, ' ') I backbone
      [ -27.06424, 744, 745, 0, 0, 1,     // 446_I_N:446_I_CA:446_I_C:446_I_O [ 446_I_N:446_I_CA:446_I_C -- 446_I_CA:446_I_C:446_I_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 154.10705, 744, 746, 0, 0, 1,     // 446_I_N:446_I_CA:446_I_C:447_S_N [ 446_I_N:446_I_CA:446_I_C -- 446_I_CA:446_I_C:447_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.02290, 746, 747, 0, 0, 1,     // 446_I_CA:446_I_C:447_S_N:447_S_CA [ 446_I_CA:446_I_C:447_S_N -- 446_I_C:447_S_N:447_S_CA ] 
        [ [ -0.411626649594191, -0.4366910850711724, -0.799915243987276, 0.0 ], [ 0.1998125671205203, -0.8996115251703716, 0.3882963324593029, 0.0 ], [ -0.8891785194012329, 0.0, 0.4575604447889167, 15.398264250410332 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.89957, 747, 748, 0, 0, 1,     // 446_I_C:447_S_N:447_S_CA:447_S_C [ 446_I_C:447_S_N:447_S_CA -- 447_S_N:447_S_CA:447_S_C ] 
        [ [ 0.8988538017416315, 0.4296082349998675, -0.08659450048885545, -10.689131561935843 ], [ -0.4259578740033765, 0.902888060982244, 0.05790543074917377, 5.188737949267897 ], [ 0.10306179053916828, -0.015162907241764405, 0.9945593765959065, 21.5125667701463 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 446, ' ') I sidechain

     ],
     [ // 187 : (' ', 447, ' ') S backbone
      [ -29.44827, 748, 749, 0, 0, 1,     // 447_S_N:447_S_CA:447_S_C:447_S_O [ 447_S_N:447_S_CA:447_S_C -- 447_S_CA:447_S_C:447_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 153.85207, 748, 750, 0, 0, 1,     // 447_S_N:447_S_CA:447_S_C:448_S_N [ 447_S_N:447_S_CA:447_S_C -- 447_S_CA:447_S_C:448_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -177.55255, 750, 751, 0, 0, 1,     // 447_S_CA:447_S_C:448_S_N:448_S_CA [ 447_S_CA:447_S_C:448_S_N -- 447_S_C:448_S_N:448_S_CA ] 
        [ [ -0.4267218439509441, -0.4406902314027435, -0.789747167035946, 0.0 ], [ 0.2094916854304306, -0.8976592448954093, 0.38771266911118746, 0.0 ], [ -0.8797850314881606, 0.0, 0.47537174755066836, 15.344778551104856 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -51.20840, 751, 752, 0, 0, 1,     // 447_S_C:448_S_N:448_S_CA:448_S_C [ 447_S_C:448_S_N:448_S_CA -- 448_S_N:448_S_CA:448_S_C ] 
        [ [ 0.904163625707469, 0.4220659232965689, -0.06594311441946396, -10.570547480707921 ], [ -0.41757179481137563, 0.905786338061457, 0.07200629110863932, 5.189426880875746 ], [ 0.0901217738703033, -0.03756948459894736, 0.9952218846576071, 21.707497960825417 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 447, ' ') S sidechain

     ],
     [ // 188 : (' ', 448, ' ') S backbone
      [ 131.60614, 752, 753, 0, 0, 1,     // 448_S_N:448_S_CA:448_S_C:448_S_O [ 448_S_N:448_S_CA:448_S_C -- 448_S_CA:448_S_C:448_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.61494, 752, 754, 0, 0, 1,     // 448_S_N:448_S_CA:448_S_C:449_Q_N [ 448_S_N:448_S_CA:448_S_C -- 448_S_CA:448_S_C:449_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.84735, 754, 755, 0, 0, 1,     // 448_S_CA:448_S_C:449_Q_N:449_Q_CA [ 448_S_CA:448_S_C:449_Q_N -- 448_S_C:449_Q_N:449_Q_CA ] 
        [ [ 0.35442773726130516, 0.6641211992249664, 0.6582735083533467, 0.0 ], [ -0.3148409919744801, 0.7476249278481774, -0.5847496190956715, 0.0 ], [ -0.8804863024671972, 0.0, 0.47407158865264587, 15.425228235445113 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.93605, 755, 756, 0, 0, 1,     // 448_S_C:449_Q_N:449_Q_CA:449_Q_C [ 448_S_C:449_Q_N:449_Q_CA -- 449_Q_N:449_Q_CA:449_Q_C ] 
        [ [ -0.7122068953931567, -0.6932684521875976, 0.1101825365284973, 8.818802889831833 ], [ 0.6977568848025154, -0.7163280165178157, 0.003082606420066392, -7.833813096335207 ], [ 0.07678976406481727, 0.07907607699591813, 0.9939065882575713, 21.77630243628169 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 448, ' ') S sidechain

     ],
     [ // 189 : (' ', 449, ' ') Q backbone
      [ 135.26361, 756, 757, 0, 0, 1,     // 449_Q_N:449_Q_CA:449_Q_C:449_Q_O [ 449_Q_N:449_Q_CA:449_Q_C -- 449_Q_CA:449_Q_C:449_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.27274, 756, 758, 0, 0, 1,     // 449_Q_N:449_Q_CA:449_Q_C:450_Q_N [ 449_Q_N:449_Q_CA:449_Q_C -- 449_Q_CA:449_Q_C:450_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.70325, 758, 759, 0, 0, 1,     // 449_Q_CA:449_Q_C:450_Q_N:450_Q_CA [ 449_Q_CA:449_Q_C:450_Q_N -- 449_Q_C:450_Q_N:450_Q_CA ] 
        [ [ 0.34431483321414047, 0.6726605806020963, 0.6549618606703518, 0.0 ], [ -0.3130031048516368, 0.7399511762974978, -0.595400128526342, 0.0 ], [ -0.8851419953781167, 0.0, 0.46532101609324084, 15.404134036034455 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -68.17713, 759, 760, 0, 0, 1,     // 449_Q_C:450_Q_N:450_Q_CA:450_Q_C [ 449_Q_C:450_Q_N:450_Q_CA -- 450_Q_N:450_Q_CA:450_Q_C ] 
        [ [ -0.7014175690060764, -0.7015733944625652, 0.12573052959357078, 8.757862608924157 ], [ 0.7062435086921167, -0.7078969338014608, -0.010101363508862045, -7.961429261900119 ], [ 0.09609110427013523, 0.08171109653384875, 0.9920130021241567, 21.626202484140144 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 449, ' ') Q sidechain

     ],
     [ // 190 : (' ', 450, ' ') Q backbone
      [ 145.47051, 760, 761, 0, 0, 1,     // 450_Q_N:450_Q_CA:450_Q_C:450_Q_O [ 450_Q_N:450_Q_CA:450_Q_C -- 450_Q_CA:450_Q_C:450_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -31.29320, 760, 762, 0, 0, 1,     // 450_Q_N:450_Q_CA:450_Q_C:451_Q_N [ 450_Q_N:450_Q_CA:450_Q_C -- 450_Q_CA:450_Q_C:451_Q_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.30935, 762, 763, 0, 0, 1,     // 450_Q_CA:450_Q_C:451_Q_N:451_Q_CA [ 450_Q_CA:450_Q_C:451_Q_N -- 450_Q_C:451_Q_N:451_Q_CA ] 
        [ [ 0.426036590244731, 0.5194177379534397, 0.7407415455285162, 0.0 ], [ -0.25896508370774035, 0.8545204582099435, -0.45025756175872406, 0.0 ], [ -0.8668505691253171, 0.0, 0.4985680402985272, 15.363410946082032 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.15282, 763, 764, 0, 0, 1,     // 450_Q_C:451_Q_N:451_Q_CA:451_Q_C [ 450_Q_C:451_Q_N:451_Q_CA -- 451_Q_N:451_Q_CA:451_Q_C ] 
        [ [ -0.8183960145423782, -0.5655172560023195, 0.10206956718218323, 9.917606599231876 ], [ 0.5721015380707121, -0.8185292524845588, 0.05205471125839491, -6.028387894277314 ], [ 0.054109089054214354, 0.10099552460716357, 0.9934143699841708, 22.03861597066699 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 450, ' ') Q sidechain

     ],
     [ // 191 : (' ', 451, ' ') Q backbone
      [ 129.47505, 764, 765, 0, 0, 1,     // 451_Q_N:451_Q_CA:451_Q_C:451_Q_O [ 451_Q_N:451_Q_CA:451_Q_C -- 451_Q_CA:451_Q_C:451_Q_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.88848, 764, 766, 0, 0, 1,     // 451_Q_N:451_Q_CA:451_Q_C:452_S_N [ 451_Q_N:451_Q_CA:451_Q_C -- 451_Q_CA:451_Q_C:452_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.04698, 766, 767, 0, 0, 1,     // 451_Q_CA:451_Q_C:452_S_N:452_S_CA [ 451_Q_CA:451_Q_C:452_S_N -- 451_Q_C:452_S_N:452_S_CA ] 
        [ [ 0.33204351681299116, 0.7179863061797072, 0.6117538451704905, 0.0 ], [ -0.3425044503001712, 0.6960572276317657, -0.6310269704110576, 0.0 ], [ -0.8788844090476507, 0.0, 0.4770347948870833, 15.352429644676338 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.69684, 767, 768, 0, 0, 1,     // 451_Q_C:452_S_N:452_S_CA:452_S_C [ 451_Q_C:452_S_N:452_S_CA -- 452_S_N:452_S_CA:452_S_C ] 
        [ [ -0.6605056281767671, -0.74397347108347, 0.10117207851392296, 8.19261741896014 ], [ 0.7478331814080149, -0.6638866016190342, 0.00033611887300123304, -8.450723424848285 ], [ 0.06691672385869904, 0.07588184575208223, 0.9948688846040346, 21.740887356704 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 451, ' ') Q sidechain

     ],
     [ // 192 : (' ', 452, ' ') S backbone
      [ 134.66418, 768, 769, 0, 0, 1,     // 452_S_N:452_S_CA:452_S_C:452_S_O [ 452_S_N:452_S_CA:452_S_C -- 452_S_CA:452_S_C:452_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.58178, 768, 770, 0, 0, 1,     // 452_S_N:452_S_CA:452_S_C:453_M_N [ 452_S_N:452_S_CA:452_S_C -- 452_S_CA:452_S_C:453_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.66343, 770, 771, 0, 0, 1,     // 452_S_CA:452_S_C:453_M_N:453_M_CA [ 452_S_CA:452_S_C:453_M_N -- 452_S_C:453_M_N:453_M_CA ] 
        [ [ 0.3464134365005736, 0.6766418593212785, 0.6497334262804294, 0.0 ], [ -0.31834023918280036, 0.7363122939447929, -0.5970793062088103, 0.0 ], [ -0.8824155614725414, 0.0, 0.47047080342046677, 15.334309465392536 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.19453, 771, 772, 0, 0, 1,     // 452_S_C:453_M_N:453_M_CA:453_M_C [ 452_S_C:453_M_N:453_M_CA -- 453_M_N:453_M_CA:453_M_C ] 
        [ [ -0.698767388677004, -0.7059275245778831, 0.1157171835332285, 8.700579713255237 ], [ 0.7106821313415824, -0.7035132196868713, -0.00024066075196421005, -7.9954884398432 ], [ 0.08157845740945535, 0.08206996874103284, 0.993282173160044, 21.6343835609375 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 452, ' ') S sidechain

     ],
     [ // 193 : (' ', 453, ' ') M backbone
      [ 141.19816, 772, 773, 0, 0, 1,     // 453_M_N:453_M_CA:453_M_C:453_M_O [ 453_M_N:453_M_CA:453_M_C -- 453_M_CA:453_M_C:453_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.95330, 772, 774, 0, 0, 1,     // 453_M_N:453_M_CA:453_M_C:454_R_N [ 453_M_N:453_M_CA:453_M_C -- 453_M_CA:453_M_C:454_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.83368, 774, 775, 0, 0, 1,     // 453_M_CA:453_M_C:454_R_N:454_R_CA [ 453_M_CA:453_M_C:454_R_N -- 453_M_C:454_R_N:454_R_CA ] 
        [ [ 0.37470651654176906, 0.6150190543301111, 0.6937914594977578, 0.0 ], [ -0.29226134215499966, 0.7885122464558785, -0.541139302833256, 0.0 ], [ -0.8798740445898442, 0.0, 0.47520697138942397, 15.343612379281536 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.24060, 775, 776, 0, 0, 1,     // 453_M_C:454_R_N:454_R_CA:454_R_C [ 453_M_C:454_R_N:454_R_CA -- 454_R_N:454_R_CA:454_R_C ] 
        [ [ -0.7652746779586111, -0.6347769443670133, 0.10683116668567115, 9.288368561758332 ], [ 0.636131858804199, -0.7711655578616331, -0.0252970469081001, -7.24468602079193 ], [ 0.0984424983919975, 0.048599519216077025, 0.9939553114915719, 21.705606942571748 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 453, ' ') M sidechain

     ],
     [ // 194 : (' ', 454, ' ') R backbone
      [ 123.85829, 776, 777, 0, 0, 1,     // 454_R_N:454_R_CA:454_R_C:454_R_O [ 454_R_N:454_R_CA:454_R_C -- 454_R_CA:454_R_C:454_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -51.78468, 776, 778, 0, 0, 1,     // 454_R_N:454_R_CA:454_R_C:455_L_N [ 454_R_N:454_R_CA:454_R_C -- 454_R_CA:454_R_C:455_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.07345, 778, 779, 0, 0, 1,     // 454_R_CA:454_R_C:455_L_N:455_L_CA [ 454_R_CA:454_R_C:455_L_N -- 454_R_C:455_L_N:455_L_CA ] 
        [ [ 0.2857109165736507, 0.7856915634987154, 0.548687560636824, 0.0 ], [ -0.36287417956864787, 0.618618434131205, -0.6968741369540389, 0.0 ], [ -0.8869563698136594, 0.0, 0.46185322132358775, 15.21527409427323 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -55.02507, 779, 780, 0, 0, 1,     // 454_R_C:455_L_N:455_L_CA:455_L_C [ 454_R_C:455_L_N:455_L_CA -- 455_L_N:455_L_CA:455_L_C ] 
        [ [ -0.5926069226690606, -0.7992540548220652, 0.10004994280454686, 7.313612727665875 ], [ 0.7984043723246852, -0.5992848192343655, -0.0583794800255988, -9.288833797676912 ], [ 0.10661844801689507, 0.04528422778099489, 0.9932682644969322, 21.371447188238175 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 454, ' ') R sidechain

     ],
     [ // 195 : (' ', 455, ' ') L backbone
      [ 126.51534, 780, 781, 0, 0, 1,     // 455_L_N:455_L_CA:455_L_C:455_L_O [ 455_L_N:455_L_CA:455_L_C -- 455_L_CA:455_L_C:455_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.82566, 780, 782, 0, 0, 1,     // 455_L_N:455_L_CA:455_L_C:456_A_N [ 455_L_N:455_L_CA:455_L_C -- 455_L_CA:455_L_C:456_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.73816, 782, 783, 0, 0, 1,     // 455_L_CA:455_L_C:456_A_N:456_A_CA [ 455_L_CA:455_L_C:456_A_N -- 455_L_C:456_A_N:456_A_CA ] 
        [ [ 0.31227823651937914, 0.7527098653076314, 0.5795775717407584, 0.0 ], [ -0.35703510131404687, 0.6583523818355693, -0.6626447598533477, 0.0 ], [ -0.8803455227500251, 0.0, 0.47433296382813733, 15.34588528861059 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.67719, 783, 784, 0, 0, 1,     // 455_L_C:456_A_N:456_A_CA:456_A_C [ 455_L_C:456_A_N:456_A_CA -- 456_A_N:456_A_CA:456_A_C ] 
        [ [ -0.6404031309440986, -0.7644478806109865, 0.0741840124714671, 7.748459411281832 ], [ 0.7644678204001277, -0.6437486289204187, -0.034302395481351534, -8.858997097491002 ], [ 0.07397824984192125, 0.034743928857458126, 0.9966544426017847, 21.68731376489969 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 455, ' ') L sidechain

     ],
     [ // 196 : (' ', 456, ' ') A backbone
      [ 131.28540, 784, 785, 0, 0, 1,     // 456_A_N:456_A_CA:456_A_C:456_A_O [ 456_A_N:456_A_CA:456_A_C -- 456_A_CA:456_A_C:456_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.56067, 784, 786, 0, 0, 1,     // 456_A_N:456_A_CA:456_A_C:457_N_N [ 456_A_N:456_A_CA:456_A_C -- 456_A_CA:456_A_C:457_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.89202, 786, 787, 0, 0, 1,     // 456_A_CA:456_A_C:457_N_N:457_N_CA [ 456_A_CA:456_A_C:457_N_N -- 456_A_C:457_N_N:457_N_CA ] 
        [ [ 0.33512841613259897, 0.7139922072716745, 0.6147390280889752, 0.0 ], [ -0.34175224100255497, 0.7001536459630288, -0.6268893664869498, 0.0 ], [ -0.8780058943254238, 0.0, 0.4786498193145098, 15.419393627425652 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.66750, 787, 788, 0, 0, 1,     // 456_A_C:457_N_N:457_N_CA:457_N_C [ 456_A_C:457_N_N:457_N_CA -- 457_N_N:457_N_CA:457_N_C ] 
        [ [ -0.6773614517170359, -0.7311119602092438, 0.08158900273190359, 8.233600165216064 ], [ 0.7324402809543312, -0.6805947396794062, -0.017945338009123343, -8.39633756054948 ], [ 0.0686490973234899, 0.04760359187832647, 0.9965044904449512, 21.830262388499627 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 456, ' ') A sidechain

     ],
     [ // 197 : (' ', 457, ' ') N backbone
      [ 141.90133, 788, 789, 0, 0, 1,     // 457_N_N:457_N_CA:457_N_C:457_N_O [ 457_N_N:457_N_CA:457_N_C -- 457_N_CA:457_N_C:457_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.81573, 788, 790, 0, 0, 1,     // 457_N_N:457_N_CA:457_N_C:458_L_N [ 457_N_N:457_N_CA:457_N_C -- 457_N_CA:457_N_C:458_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.35173, 790, 791, 0, 0, 1,     // 457_N_CA:457_N_C:458_L_N:458_L_CA [ 457_N_CA:457_N_C:458_L_N -- 457_N_C:458_L_N:458_L_CA ] 
        [ [ 0.3765052529480942, 0.6131239079334398, 0.6944946853813343, 0.0 ], [ -0.2922129649663589, 0.7899867552816483, -0.5390106766893794, 0.0 ], [ -0.8791219355743896, 0.0, 0.4765969181519525, 15.37169956862465 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.50965, 791, 792, 0, 0, 1,     // 457_N_C:458_L_N:458_L_CA:458_L_C [ 457_N_C:458_L_N:458_L_CA -- 458_L_N:458_L_CA:458_L_C ] 
        [ [ -0.777639523484059, -0.6236999572693221, 0.07921574854617966, 9.285169608948262 ], [ 0.6234425761930199, -0.7812547181582643, -0.030990636412661702, -7.206398636941069 ], [ 0.08121653591046772, 0.025286926616278748, 0.9963756548797279, 21.743646414716164 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 457, ' ') N sidechain

     ],
     [ // 198 : (' ', 458, ' ') L backbone
      [ 135.49044, 792, 793, 0, 0, 1,     // 458_L_N:458_L_CA:458_L_C:458_L_O [ 458_L_N:458_L_CA:458_L_C -- 458_L_CA:458_L_C:458_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.62092, 792, 794, 0, 0, 1,     // 458_L_N:458_L_CA:458_L_C:459_L_N [ 458_L_N:458_L_CA:458_L_C -- 458_L_CA:458_L_C:459_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.83250, 794, 795, 0, 0, 1,     // 458_L_CA:458_L_C:459_L_N:459_L_CA [ 458_L_CA:458_L_C:459_L_N -- 458_L_C:459_L_N:459_L_CA ] 
        [ [ 0.36109338934552737, 0.6510514214625815, 0.6676403304044113, 0.0 ], [ -0.30972325350800306, 0.7590336267989397, -0.5726599860518152, 0.0 ], [ -0.879592559317879, 0.0, 0.4757277893844581, 15.328563528375705 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.65347, 795, 796, 0, 0, 1,     // 458_L_C:459_L_N:459_L_CA:459_L_C [ 458_L_C:459_L_N:459_L_CA -- 459_L_N:459_L_CA:459_L_C ] 
        [ [ -0.7308206027302074, -0.6755715040316481, 0.09749045884327918, 8.920686977455288 ], [ 0.6785820253920789, -0.7345182762363751, -0.0030555997002371387, -7.651605583784551 ], [ 0.07367279986426656, 0.06392217780364595, 0.9952317688583895, 21.685008012496233 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 458, ' ') L sidechain

     ],
     [ // 199 : (' ', 459, ' ') L backbone
      [ 143.64993, 796, 797, 0, 0, 1,     // 459_L_N:459_L_CA:459_L_C:459_L_O [ 459_L_N:459_L_CA:459_L_C -- 459_L_CA:459_L_C:459_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -31.86779, 796, 798, 0, 0, 1,     // 459_L_N:459_L_CA:459_L_C:460_M_N [ 459_L_N:459_L_CA:459_L_C -- 459_L_CA:459_L_C:460_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.33238, 798, 799, 0, 0, 1,     // 459_L_CA:459_L_C:460_M_N:460_M_CA [ 459_L_CA:459_L_C:460_M_N -- 459_L_C:460_M_N:460_M_CA ] 
        [ [ 0.4115362388260554, 0.527960996199827, 0.7428964333099132, 0.0 ], [ -0.2558378796430702, 0.8492686185722903, -0.46183307898855797, 0.0 ], [ -0.874748480120224, 0.0, 0.48457723484224696, 15.389608153230315 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.04100, 799, 800, 0, 0, 1,     // 459_L_C:460_M_N:460_M_CA:460_M_C [ 459_L_C:460_M_N:460_M_CA -- 460_M_N:460_M_CA:460_M_C ] 
        [ [ -0.8072015552569478, -0.5781500382330393, 0.11903017466970252, 9.941537612379854 ], [ 0.5869506996479478, -0.8075396687297602, 0.05803929367776303, -6.180310901412803 ], [ 0.06256616796279978, 0.11671425242430436, 0.9911928459727131, 21.8742843345703 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 459, ' ') L sidechain

     ],
     [ // 200 : (' ', 460, ' ') M backbone
      [ 155.41755, 800, 801, 0, 0, 1,     // 460_M_N:460_M_CA:460_M_C:460_M_O [ 460_M_N:460_M_CA:460_M_C -- 460_M_CA:460_M_C:460_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -21.97238, 800, 802, 0, 0, 1,     // 460_M_N:460_M_CA:460_M_C:461_L_N [ 460_M_N:460_M_CA:460_M_C -- 460_M_CA:460_M_C:461_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.30427, 802, 803, 0, 0, 1,     // 460_M_CA:460_M_C:461_L_N:461_L_CA [ 460_M_CA:460_M_C:461_L_N -- 460_M_C:461_L_N:461_L_CA ] 
        [ [ 0.4571222691783536, 0.37415965397128936, 0.80687284274619, 0.0 ], [ -0.1844331394043388, 0.9273643045427643, -0.3255454864521619, 0.0 ], [ -0.8700710592306198, 0.0, 0.4929263148679602, 15.400582255322897 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.81217, 803, 804, 0, 0, 1,     // 460_M_C:461_L_N:461_L_CA:461_L_C [ 460_M_C:461_L_N:461_L_CA -- 461_L_N:461_L_CA:461_L_C ] 
        [ [ -0.9060867211729694, -0.41767976230198695, 0.06745717068906164, 10.811105989847364 ], [ 0.42169761490391516, -0.904481718116808, 0.06390573665029052, -4.361910045915176 ], [ 0.03432164475031453, 0.08635066737334858, 0.9956734338857347, 22.00519001121772 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 460, ' ') M sidechain

     ],
     [ // 201 : (' ', 461, ' ') L backbone
      [ 142.75865, 804, 805, 0, 0, 1,     // 461_L_N:461_L_CA:461_L_C:461_L_O [ 461_L_N:461_L_CA:461_L_C -- 461_L_CA:461_L_C:461_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -32.42730, 804, 806, 0, 0, 1,     // 461_L_N:461_L_CA:461_L_C:462_L_N [ 461_L_N:461_L_CA:461_L_C -- 461_L_CA:461_L_C:462_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.70018, 806, 807, 0, 0, 1,     // 461_L_CA:461_L_C:462_L_N:462_L_CA [ 461_L_CA:461_L_C:462_L_N -- 461_L_C:462_L_N:462_L_CA ] 
        [ [ 0.40105236273855216, 0.5362290784942757, 0.742708137641567, 0.0 ], [ -0.25478373053291875, 0.844072494147973, -0.4718335249625117, 0.0 ], [ -0.8799103664564669, 0.0, 0.4751397131396682, 15.347644437069013 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.29461, 807, 808, 0, 0, 1,     // 461_L_C:462_L_N:462_L_CA:462_L_C [ 461_L_C:462_L_N:462_L_CA -- 462_L_N:462_L_CA:462_L_C ] 
        [ [ -0.8281915409814697, -0.5518908523731961, 0.09754618656577724, 9.95368009170841 ], [ 0.5525634735078861, -0.8331684807506313, -0.022447503838609856, -6.323452950082647 ], [ 0.09366098009116726, 0.03530962688094842, 0.9949778143546167, 21.715406501015558 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 461, ' ') L sidechain

     ],
     [ // 202 : (' ', 462, ' ') L backbone
      [ 138.08358, 808, 809, 0, 0, 1,     // 462_L_N:462_L_CA:462_L_C:462_L_O [ 462_L_N:462_L_CA:462_L_C -- 462_L_CA:462_L_C:462_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.93887, 808, 810, 0, 0, 1,     // 462_L_N:462_L_CA:462_L_C:463_S_N [ 462_L_N:462_L_CA:462_L_C -- 462_L_CA:462_L_C:463_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.90305, 810, 811, 0, 0, 1,     // 462_L_CA:462_L_C:463_S_N:463_S_CA [ 462_L_CA:462_L_C:463_S_N -- 462_L_C:463_S_N:463_S_CA ] 
        [ [ 0.3772729056759267, 0.6284909024462568, 0.6801943400125681, 0.0 ], [ -0.30484369529394356, 0.7778169357517807, -0.5496099852597691, 0.0 ], [ -0.8744915528936666, 0.0, 0.48504074459536217, 15.38167427836341 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -70.34235, 811, 812, 0, 0, 1,     // 462_L_C:463_S_N:463_S_CA:463_S_C [ 462_L_C:463_S_N:463_S_CA -- 463_S_N:463_S_CA:463_S_C ] 
        [ [ -0.7647937275432639, -0.6418746089184894, 0.05556564348782704, 9.123923967877385 ], [ 0.6426533567121192, -0.766141705564434, -0.004852844677702898, -7.372304387895507 ], [ 0.045686074652191815, 0.031998022135072066, 0.9984432428347231, 21.887866133663103 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 462, ' ') L sidechain

     ],
     [ // 203 : (' ', 463, ' ') S backbone
      [ 142.22629, 812, 813, 0, 0, 1,     // 463_S_N:463_S_CA:463_S_C:463_S_O [ 463_S_N:463_S_CA:463_S_C -- 463_S_CA:463_S_C:463_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.41006, 812, 814, 0, 0, 1,     // 463_S_N:463_S_CA:463_S_C:464_H_N [ 463_S_N:463_S_CA:463_S_C -- 463_S_CA:463_S_C:464_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.03611, 814, 815, 0, 0, 1,     // 463_S_CA:463_S_C:464_H_N:464_H_CA [ 463_S_CA:463_S_C:464_H_N -- 463_S_C:464_H_N:464_H_CA ] 
        [ [ 0.38757731656217637, 0.5794243017702044, 0.7169737109577824, 0.0 ], [ -0.2755393073252384, 0.815026060023918, -0.5097162069236338, 0.0 ], [ -0.8796942160803323, 0.0, 0.4755397840294854, 15.299130225889737 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -70.13488, 815, 816, 0, 0, 1,     // 463_S_C:464_H_N:464_H_CA:464_H_C [ 463_S_C:464_H_N:464_H_CA -- 464_H_N:464_H_CA:464_H_C ] 
        [ [ -0.7911867795502371, -0.6048305426406372, 0.09057314477237154, 9.591098531471992 ], [ 0.6078363554529032, -0.7940289717790705, 0.0072771536471400765, -6.81857408294933 ], [ 0.06751625622512406, 0.060811237978719505, 0.9958632177567555, 21.660519569118065 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 463, ' ') S sidechain

     ],
     [ // 204 : (' ', 464, ' ') H backbone
      [ 131.49716, 816, 817, 0, 0, 1,     // 464_H_N:464_H_CA:464_H_C:464_H_O [ 464_H_N:464_H_CA:464_H_C -- 464_H_CA:464_H_C:464_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.05556, 816, 818, 0, 0, 1,     // 464_H_N:464_H_CA:464_H_C:465_V_N [ 464_H_N:464_H_CA:464_H_C -- 464_H_CA:464_H_C:465_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.99220, 818, 819, 0, 0, 1,     // 464_H_CA:464_H_C:465_V_N:465_V_CA [ 464_H_CA:464_H_C:465_V_N -- 464_H_C:465_V_N:465_V_CA ] 
        [ [ 0.33271501919503976, 0.6953556058100027, 0.6370096525725873, 0.0 ], [ -0.3219232676369267, 0.7186658343546075, -0.6163479766214532, 0.0 ], [ -0.8863780941313971, 0.0, 0.4629620656641233, 15.349756236424135 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.11278, 819, 820, 0, 0, 1,     // 464_H_C:465_V_N:465_V_CA:465_V_C [ 464_H_C:465_V_N:465_V_CA -- 465_V_N:465_V_CA:465_V_C ] 
        [ [ -0.6798376705952067, -0.721744456975007, 0.13002184610878162, 8.538982839263502 ], [ 0.7254515398836434, -0.6878214411037118, -0.024934482919516377, -8.262017340130734 ], [ 0.10742813840020614, 0.07737314769262975, 0.9911975540203883, 21.55566711700606 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 464, ' ') H sidechain

     ],
     [ // 205 : (' ', 465, ' ') V backbone
      [ 127.29273, 820, 821, 0, 0, 1,     // 465_V_N:465_V_CA:465_V_C:465_V_O [ 465_V_N:465_V_CA:465_V_C -- 465_V_CA:465_V_C:465_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.31843, 820, 822, 0, 0, 1,     // 465_V_N:465_V_CA:465_V_C:466_R_N [ 465_V_N:465_V_CA:465_V_C -- 465_V_CA:465_V_C:466_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.12783, 822, 823, 0, 0, 1,     // 465_V_CA:465_V_C:466_R_N:466_R_CA [ 465_V_CA:465_V_C:466_R_N -- 465_V_C:466_R_N:466_R_CA ] 
        [ [ 0.31859260079114915, 0.7468521539954617, 0.5837042185846105, 0.0 ], [ -0.35781218831803335, 0.6649901202817522, -0.655559743897321, 0.0 ], [ -0.8777637453279736, 0.0, 0.4790937354921355, 15.325555322190947 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.14386, 823, 824, 0, 0, 1,     // 465_V_C:466_R_N:466_R_CA:466_R_C [ 465_V_C:466_R_N:466_R_CA -- 466_R_N:466_R_CA:466_R_C ] 
        [ [ -0.6340579048010588, -0.7666620140220404, 0.10099470092514874, 7.837156405499129 ], [ 0.7684065620483461, -0.6393087782154975, -0.028907464393284708, -8.801930982322551 ], [ 0.08672905372672249, 0.05927598461667716, 0.9944669068839801, 21.758149808232734 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 465, ' ') V sidechain

     ],
     [ // 206 : (' ', 466, ' ') R backbone
      [ 136.33804, 824, 825, 0, 0, 1,     // 466_R_N:466_R_CA:466_R_C:466_R_O [ 466_R_N:466_R_CA:466_R_C -- 466_R_CA:466_R_C:466_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.66993, 824, 826, 0, 0, 1,     // 466_R_N:466_R_CA:466_R_C:467_H_N [ 466_R_N:466_R_CA:466_R_C -- 466_R_CA:466_R_C:467_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.26560, 826, 827, 0, 0, 1,     // 466_R_CA:466_R_C:467_H_N:467_H_CA [ 466_R_CA:466_R_C:467_H_N -- 466_R_C:467_H_N:467_H_CA ] 
        [ [ 0.33914114127178174, 0.6905029499596353, 0.6388966758357056, 0.0 ], [ -0.32375001125474046, 0.7233295764014087, -0.609901839737768, 0.0 ], [ -0.8832718814212466, 0.0, 0.4688611558772716, 15.277267369570453 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -66.59815, 827, 828, 0, 0, 1,     // 466_R_C:467_H_N:467_H_CA:467_H_C [ 466_R_C:467_H_N:467_H_CA -- 467_H_N:467_H_CA:467_H_C ] 
        [ [ -0.7259763101992908, -0.6860993084099472, 0.0471819460049643, 8.513122791903708 ], [ 0.6835627489515372, -0.7274198015245538, -0.060020001631093875, -8.126774561636283 ], [ 0.07550076340833074, -0.011321278610260544, 0.997081472787148, 21.52471354248424 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 466, ' ') R sidechain

     ],
     [ // 207 : (' ', 467, ' ') H backbone
      [ 129.99012, 828, 829, 0, 0, 1,     // 467_H_N:467_H_CA:467_H_C:467_H_O [ 467_H_N:467_H_CA:467_H_C -- 467_H_CA:467_H_C:467_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.65633, 828, 830, 0, 0, 1,     // 467_H_N:467_H_CA:467_H_C:468_A_N [ 467_H_N:467_H_CA:467_H_C -- 467_H_CA:467_H_C:468_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.91638, 830, 831, 0, 0, 1,     // 467_H_CA:467_H_C:468_A_N:468_A_CA [ 467_H_CA:467_H_C:468_A_N -- 467_H_C:468_A_N:468_A_CA ] 
        [ [ 0.3233840238534774, 0.7272498273088985, 0.6054175929022139, 0.0 ], [ -0.3426431794849863, 0.6863728496081246, -0.6414732752603511, 0.0 ], [ -0.882053527099548, 0.0, 0.47114920707908126, 15.365914800406529 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.68904, 831, 832, 0, 0, 1,     // 467_H_C:468_A_N:468_A_CA:468_A_C [ 467_H_C:468_A_N:468_A_CA -- 468_A_N:468_A_CA:468_A_C ] 
        [ [ -0.6473718393649426, -0.7530440376345334, 0.11761963688232654, 8.102399186759532 ], [ 0.7569030259643109, -0.6533115319061169, -0.016788435440871066, -8.584938073705649 ], [ 0.0894846963637942, 0.07815829873764014, 0.9929167988583522, 21.67137899080136 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 467, ' ') H sidechain

     ],
     [ // 208 : (' ', 468, ' ') A backbone
      [ 133.53215, 832, 833, 0, 0, 1,     // 468_A_N:468_A_CA:468_A_C:468_A_O [ 468_A_N:468_A_CA:468_A_C -- 468_A_CA:468_A_C:468_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -41.77153, 832, 834, 0, 0, 1,     // 468_A_N:468_A_CA:468_A_C:469_S_N [ 468_A_N:468_A_CA:468_A_C -- 468_A_CA:468_A_C:469_S_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.76754, 834, 835, 0, 0, 1,     // 468_A_CA:468_A_C:469_S_N:469_S_CA [ 468_A_CA:468_A_C:469_S_N -- 468_A_C:469_S_N:469_S_CA ] 
        [ [ 0.3446533186578109, 0.6661619847176554, 0.6613942092694732, 0.0 ], [ -0.3078476230447169, 0.7458070863950237, -0.5907635998169518, 0.0 ], [ -0.8868167403268137, 0.0, 0.46212127096263883, 15.376973403921571 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.32400, 835, 836, 0, 0, 1,     // 468_A_C:469_S_N:469_S_CA:469_S_C [ 468_A_C:469_S_N:469_S_CA -- 469_S_N:469_S_CA:469_S_C ] 
        [ [ -0.7192137626792194, -0.6845361456751122, 0.1189194216142208, 8.838823772089334 ], [ 0.6849656998687615, -0.7272618184114467, -0.04372913770214785, -7.894921480359793 ], [ 0.11641973018392024, 0.05000512718849361, 0.9919404889804457, 21.552728223251275 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 468, ' ') A sidechain

     ],
     [ // 209 : (' ', 469, ' ') S backbone
      [ 133.88384, 836, 837, 0, 0, 1,     // 469_S_N:469_S_CA:469_S_C:469_S_O [ 469_S_N:469_S_CA:469_S_C -- 469_S_CA:469_S_C:469_S_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -44.01826, 836, 838, 0, 0, 1,     // 469_S_N:469_S_CA:469_S_C:470_N_N [ 469_S_N:469_S_CA:469_S_C -- 469_S_CA:469_S_C:470_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.56545, 838, 839, 0, 0, 1,     // 469_S_CA:469_S_C:470_N_N:470_N_CA [ 469_S_CA:469_S_C:470_N_N -- 469_S_C:470_N_N:470_N_CA ] 
        [ [ 0.34601395426044296, 0.6948875281497036, 0.630401195017147, 0.0 ], [ -0.3343549125662301, 0.7191184347671772, -0.6091596418188204, 0.0 ], [ -0.8766305583881278, 0.0, 0.4811640719132293, 15.273696882879205 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.90804, 839, 840, 0, 0, 1,     // 469_S_C:470_N_N:470_N_CA:470_N_C [ 469_S_C:470_N_N:470_N_CA -- 470_N_N:470_N_CA:470_N_C ] 
        [ [ -0.695366128831881, -0.7143685458607749, 0.0783806580613231, 8.43009788639625 ], [ 0.716277411569873, -0.6977961708688449, -0.005212829897541974, -8.146043265090979 ], [ 0.05841760477909815, 0.05251746952719496, 0.9969098749868682, 21.708107632267453 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 469, ' ') S sidechain

     ],
     [ // 210 : (' ', 470, ' ') N backbone
      [ 128.48811, 840, 841, 0, 0, 1,     // 470_N_N:470_N_CA:470_N_C:470_N_O [ 470_N_N:470_N_CA:470_N_C -- 470_N_CA:470_N_C:470_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.37616, 840, 842, 0, 0, 1,     // 470_N_N:470_N_CA:470_N_C:471_K_N [ 470_N_N:470_N_CA:470_N_C -- 470_N_CA:470_N_C:471_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.55851, 842, 843, 0, 0, 1,     // 470_N_CA:470_N_C:471_K_N:471_K_CA [ 470_N_CA:470_N_C:471_K_N -- 470_N_C:471_K_N:471_K_CA ] 
        [ [ 0.3236424286763155, 0.7238848943778956, 0.6092997932479536, 0.0 ], [ -0.3395750320129967, 0.6899207633427935, -0.639295032001555, 0.0 ], [ -0.8831445951789935, 0.0, 0.46910086762457787, 15.394704403506484 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.26217, 843, 844, 0, 0, 1,     // 470_N_C:471_K_N:471_K_CA:471_K_C [ 470_N_C:471_K_N:471_K_CA -- 471_K_N:471_K_CA:471_K_C ] 
        [ [ -0.6554990981832222, -0.7467742370683517, 0.11246853396378872, 8.167479397849775 ], [ 0.7495068358897052, -0.6615518256708757, -0.024262829722309025, -8.569556498922955 ], [ 0.09252262012926796, 0.06839167202600649, 0.993359020677772, 21.68285961756549 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 470, ' ') N sidechain

     ],
     [ // 211 : (' ', 471, ' ') K backbone
      [ 134.54154, 844, 845, 0, 0, 1,     // 471_K_N:471_K_CA:471_K_C:471_K_O [ 471_K_N:471_K_CA:471_K_C -- 471_K_CA:471_K_C:471_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.14106, 844, 846, 0, 0, 1,     // 471_K_N:471_K_CA:471_K_C:472_G_N [ 471_K_N:471_K_CA:471_K_C -- 471_K_CA:471_K_C:472_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 173.63778, 846, 847, 0, 0, 1,     // 471_K_CA:471_K_C:472_G_N:472_G_CA [ 471_K_CA:471_K_C:472_G_N -- 471_K_C:472_G_N:472_G_CA ] 
        [ [ 0.34757049710905796, 0.6709581990443522, 0.6549884309470831, 0.0 ], [ -0.31450679993538827, 0.7414951753957403, -0.5926806709006479, 0.0 ], [ -0.8833347170431851, 0.0, 0.468742762788116, 15.349960187034231 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.71895, 847, 848, 0, 0, 1,     // 471_K_C:472_G_N:472_G_CA:472_G_C [ 471_K_C:472_G_N:472_G_CA -- 472_G_N:472_G_CA:472_G_C ] 
        [ [ -0.6955234419272845, -0.7053414605618006, 0.1368961860030044, 8.745959817644206 ], [ 0.7121009925424396, -0.702076807018652, 0.0005774656413725473, -7.913973877212282 ], [ 0.09570432670311084, 0.09788555081849604, 0.9905852314623229, 21.609010231977088 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 471, ' ') K sidechain

     ],
     [ // 212 : (' ', 472, ' ') G backbone
      [ 137.67526, 848, 849, 0, 0, 1,     // 472_G_N:472_G_CA:472_G_C:472_G_O [ 472_G_N:472_G_CA:472_G_C -- 472_G_CA:472_G_C:472_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -39.64194, 848, 850, 0, 0, 1,     // 472_G_N:472_G_CA:472_G_C:473_M_N [ 472_G_N:472_G_CA:472_G_C -- 472_G_CA:472_G_C:473_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.66290, 850, 851, 0, 0, 1,     // 472_G_CA:472_G_C:473_M_N:473_M_CA [ 472_G_CA:472_G_C:473_M_N -- 472_G_C:473_M_N:473_M_CA ] 
        [ [ 0.3559670968141077, 0.6379878536963163, 0.6828315491552097, 0.0 ], [ -0.2949207688650699, 0.7700464262211518, -0.5657298308875609, 0.0 ], [ -0.8867407547179569, 0.0, 0.4622670590927157, 15.311217239105746 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.99584, 851, 852, 0, 0, 1,     // 472_G_C:473_M_N:473_M_CA:473_M_C [ 472_G_C:473_M_N:473_M_CA -- 473_M_N:473_M_CA:473_M_C ] 
        [ [ -0.7370832082559731, -0.6630807124911557, 0.13050790332122839, 9.12653618050568 ], [ 0.6659187999930872, -0.745538074243489, -0.02692826895804426, -7.561387250449209 ], [ 0.11515422668251749, 0.06705929149291123, 0.9910815080010423, 21.4897499168726 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 472, ' ') G sidechain

     ],
     [ // 213 : (' ', 473, ' ') M backbone
      [ 128.99309, 852, 853, 0, 0, 1,     // 473_M_N:473_M_CA:473_M_C:473_M_O [ 473_M_N:473_M_CA:473_M_C -- 473_M_CA:473_M_C:473_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -49.02798, 852, 854, 0, 0, 1,     // 473_M_N:473_M_CA:473_M_C:474_E_N [ 473_M_N:473_M_CA:473_M_C -- 473_M_CA:473_M_C:474_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.18683, 854, 855, 0, 0, 1,     // 473_M_CA:473_M_C:474_E_N:474_E_CA [ 473_M_CA:473_M_C:474_E_N -- 473_M_C:474_E_N:474_E_CA ] 
        [ [ 0.30178545801579687, 0.7550298828257053, 0.582112887136505, 0.0 ], [ -0.3475070635869793, 0.6556903812320276, -0.6703051280700106, 0.0 ], [ -0.8877862231907808, 0.0, 0.46025603951784183, 15.321418603824387 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.76362, 855, 856, 0, 0, 1,     // 473_M_C:474_E_N:474_E_CA:474_E_C [ 473_M_C:474_E_N:474_E_CA -- 474_E_N:474_E_CA:474_E_C ] 
        [ [ -0.6290025565942116, -0.7689314210863321, 0.11445633868030076, 7.781550731691436 ], [ 0.7669510080455896, -0.6378447710136576, -0.07028655168941475, -8.960484254949266 ], [ 0.12705091521039769, 0.04357198362123101, 0.9909387201979301, 21.474015119413302 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 473, ' ') M sidechain

     ],
     [ // 214 : (' ', 474, ' ') E backbone
      [ 129.82168, 856, 857, 0, 0, 1,     // 474_E_N:474_E_CA:474_E_C:474_E_O [ 474_E_N:474_E_CA:474_E_C -- 474_E_CA:474_E_C:474_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.29051, 856, 858, 0, 0, 1,     // 474_E_N:474_E_CA:474_E_C:475_H_N [ 474_E_N:474_E_CA:474_E_C -- 474_E_CA:474_E_C:475_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.66505, 858, 859, 0, 0, 1,     // 474_E_CA:474_E_C:475_H_N:475_H_CA [ 474_E_CA:474_E_C:475_H_N -- 474_E_C:475_H_N:475_H_CA ] 
        [ [ 0.32266175241087197, 0.7228526583619096, 0.6110429017919, 0.0 ], [ -0.33753424659594905, 0.6910021955819825, -0.6392077894360993, 0.0 ], [ -0.8842850365725128, 0.0, 0.46694750678630864, 15.343518788721973 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.52404, 859, 860, 0, 0, 1,     // 474_E_C:475_H_N:475_H_CA:475_H_C [ 474_E_C:475_H_N:475_H_CA -- 475_H_N:475_H_CA:475_H_C ] 
        [ [ -0.6652483190867793, -0.7403987290708381, 0.09619977101065401, 8.157824238290738 ], [ 0.7415806308122611, -0.6701966057425878, -0.029911162719013062, -8.533843994708485 ], [ 0.08661894686674282, 0.05144153614932129, 0.99491252198487, 21.577574822594507 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 474, ' ') E sidechain

     ],
     [ // 215 : (' ', 475, ' ') H backbone
      [ 130.44938, 860, 861, 0, 0, 1,     // 475_H_N:475_H_CA:475_H_C:475_H_O [ 475_H_N:475_H_CA:475_H_C -- 475_H_CA:475_H_C:475_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.90665, 860, 862, 0, 0, 1,     // 475_H_N:475_H_CA:475_H_C:476_L_N [ 475_H_N:475_H_CA:475_H_C -- 475_H_CA:475_H_C:476_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.47664, 862, 863, 0, 0, 1,     // 475_H_CA:475_H_C:476_L_N:476_L_CA [ 475_H_CA:475_H_C:476_L_N -- 475_H_C:476_L_N:476_L_CA ] 
        [ [ 0.336867855136715, 0.7182071237966335, 0.6088502077714653, 0.0 ], [ -0.34770146242363026, 0.6958293809029389, -0.6284307167068631, 0.0 ], [ -0.8749992806877375, 0.0, 0.48412421835304, 15.392016008237295 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.58938, 863, 864, 0, 0, 1,     // 475_H_C:476_L_N:476_L_CA:476_L_C [ 475_H_C:476_L_N:476_L_CA -- 476_L_N:476_L_CA:476_L_C ] 
        [ [ -0.6843741013953314, -0.7269087788607093, 0.056883359206055466, 8.149990737106846 ], [ 0.7268008234769664, -0.6863399522687847, -0.02642031250545345, -8.412092916615718 ], [ 0.05824647914282974, 0.023261494683593644, 0.9980311871542638, 21.872440659200656 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 475, ' ') H sidechain

     ],
     [ // 216 : (' ', 476, ' ') L backbone
      [ 139.20137, 864, 865, 0, 0, 1,     // 476_L_N:476_L_CA:476_L_C:476_L_O [ 476_L_N:476_L_CA:476_L_C -- 476_L_CA:476_L_C:476_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.64569, 864, 866, 0, 0, 1,     // 476_L_N:476_L_CA:476_L_C:477_L_N [ 476_L_N:476_L_CA:476_L_C -- 476_L_CA:476_L_C:477_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 171.21608, 866, 867, 0, 0, 1,     // 476_L_CA:476_L_C:477_L_N:477_L_CA [ 476_L_CA:476_L_C:477_L_N -- 476_L_C:477_L_N:477_L_CA ] 
        [ [ 0.376464930345379, 0.6107767965391018, 0.696581553753386, 0.0 ], [ -0.29039558654258063, 0.7918028193997749, -0.5373255051708954, 0.0 ], [ -0.8797411889508411, 0.0, 0.4754528793301819, 15.411089398335964 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.94079, 867, 868, 0, 0, 1,     // 476_L_C:477_L_N:477_L_CA:477_L_C [ 476_L_C:477_L_N:477_L_CA -- 477_L_N:477_L_CA:477_L_C ] 
        [ [ -0.7346713410449477, -0.6611025774596765, 0.15232006671279483, 9.323223233076787 ], [ 0.6732247484276456, -0.7381700802586338, 0.04328245274352319, -7.191700105380032 ], [ 0.0838239748028644, 0.13434401619400227, 0.9873830191780262, 21.774670671654828 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 476, ' ') L sidechain

     ],
     [ // 217 : (' ', 477, ' ') L backbone
      [ 130.05286, 868, 869, 0, 0, 1,     // 477_L_N:477_L_CA:477_L_C:477_L_O [ 477_L_N:477_L_CA:477_L_C -- 477_L_CA:477_L_C:477_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.27728, 868, 870, 0, 0, 1,     // 477_L_N:477_L_CA:477_L_C:478_N_N [ 477_L_N:477_L_CA:477_L_C -- 477_L_CA:477_L_C:478_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.86915, 870, 871, 0, 0, 1,     // 477_L_CA:477_L_C:478_N_N:478_N_CA [ 477_L_CA:477_L_C:478_N_N -- 477_L_C:478_N_N:478_N_CA ] 
        [ [ 0.3260945479999718, 0.7226931670973032, 0.6094070331031333, 0.0 ], [ -0.34096769639363644, 0.6911689997611794, -0.6372020431426344, 0.0 ], [ -0.8817048121569435, 0.0, 0.47180146695330333, 15.325320974079315 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.69072, 871, 872, 0, 0, 1,     // 477_L_C:478_N_N:478_N_CA:478_N_C [ 477_L_C:478_N_N:478_N_CA -- 478_N_N:478_N_CA:478_N_C ] 
        [ [ -0.6740063937598472, -0.7343181968390479, 0.08057398439878233, 8.128320247054729 ], [ 0.7340085588088838, -0.6780133042708891, -0.0391074778907561, -8.49905233677272 ], [ 0.08334756604815244, 0.032783304023848235, 0.9959811435017878, 21.618247061425024 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 477, ' ') L sidechain

     ],
     [ // 218 : (' ', 478, ' ') N backbone
      [ 133.34333, 872, 873, 0, 0, 1,     // 478_N_N:478_N_CA:478_N_C:478_N_O [ 478_N_N:478_N_CA:478_N_C -- 478_N_CA:478_N_C:478_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.32120, 872, 874, 0, 0, 1,     // 478_N_N:478_N_CA:478_N_C:479_M_N [ 478_N_N:478_N_CA:478_N_C -- 478_N_CA:478_N_C:479_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.01909, 874, 875, 0, 0, 1,     // 478_N_CA:478_N_C:479_M_N:479_M_CA [ 478_N_CA:478_N_C:479_M_N -- 478_N_C:479_M_N:479_M_CA ] 
        [ [ 0.3479326039691419, 0.6860875825055466, 0.6389262337914661, 0.0 ], [ -0.32811824030311715, 0.7275189544827646, -0.6025401158833252, 0.0 ], [ -0.8782262370685803, 0.0, 0.4782454145356355, 15.385662727897753 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.88867, 875, 876, 0, 0, 1,     // 478_N_C:479_M_N:479_M_CA:479_M_C [ 478_N_C:479_M_N:479_M_CA -- 479_M_N:479_M_CA:479_M_C ] 
        [ [ -0.7053004148761088, -0.7032528661390727, 0.08936851259122704, 8.540495527991292 ], [ 0.7043086424451478, -0.7094713545710075, -0.02448945120358159, -8.0541240803339 ], [ 0.0806266764331857, 0.0456705956864949, 0.9956975121677161, 21.778345651380633 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 478, ' ') N sidechain

     ],
     [ // 219 : (' ', 479, ' ') M backbone
      [ 138.12456, 876, 877, 0, 0, 1,     // 479_M_N:479_M_CA:479_M_C:479_M_O [ 479_M_N:479_M_CA:479_M_C -- 479_M_CA:479_M_C:479_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -37.78072, 876, 878, 0, 0, 1,     // 479_M_N:479_M_CA:479_M_C:480_K_N [ 479_M_N:479_M_CA:479_M_C -- 479_M_CA:479_M_C:480_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 170.53427, 878, 879, 0, 0, 1,     // 479_M_CA:479_M_C:480_K_N:480_K_CA [ 479_M_CA:479_M_C:480_K_N -- 479_M_C:480_K_N:480_K_CA ] 
        [ [ 0.3676205694585912, 0.6126410928080782, 0.6996613525941424, 0.0 ], [ -0.28495763165692645, 0.7903612410808893, -0.5423359261173389, 0.0 ], [ -0.8852424894182478, 0.0, 0.4651298043864567, 15.349277176282717 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.85423, 879, 880, 0, 0, 1,     // 479_M_C:480_K_N:480_K_CA:480_K_C [ 479_M_C:480_K_N:480_K_CA -- 480_K_N:480_K_CA:480_K_C ] 
        [ [ -0.7274123178914, -0.6647574660277816, 0.17020232413273134, 9.385383684408122 ], [ 0.6792640461207017, -0.7327363492548088, 0.04120434599005323, -7.274992013747331 ], [ 0.09732253299002144, 0.1455848681733701, 0.9845472922776989, 21.588612336120928 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 479, ' ') M sidechain

     ],
     [ // 220 : (' ', 480, ' ') K backbone
      [ 129.75619, 880, 881, 0, 0, 1,     // 480_K_N:480_K_CA:480_K_C:480_K_O [ 480_K_N:480_K_CA:480_K_C -- 480_K_CA:480_K_C:480_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.81556, 880, 882, 0, 0, 1,     // 480_K_N:480_K_CA:480_K_C:481_C_N [ 480_K_N:480_K_CA:480_K_C -- 480_K_CA:480_K_C:481_C_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 179.51853, 882, 883, 0, 0, 1,     // 480_K_CA:480_K_C:481_C_N:481_C_CA [ 480_K_CA:480_K_C:481_C_N -- 480_K_C:481_C_N:481_C_CA ] 
        [ [ 0.31996132844583236, 0.7525937076951569, 0.5755236393379816, 0.0 ], [ -0.3656891539528733, 0.658485163946506, -0.6577755936058854, 0.0 ], [ -0.8740115508277965, 0.0, 0.4859051440554937, 15.335270500796149 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -62.76446, 883, 884, 0, 0, 1,     // 480_K_C:481_C_N:481_C_CA:481_C_C [ 480_K_C:481_C_N:481_C_CA -- 481_C_N:481_C_CA:481_C_C ] 
        [ [ -0.6532631864993945, -0.7552558405037021, 0.05325246050375462, 7.681294068407525 ], [ 0.75362206894525, -0.6553889482971635, -0.05019067292454714, -8.779079468082386 ], [ 0.07280787294886723, 0.0073445105340262185, 0.9973189418644756, 21.82046076185212 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 480, ' ') K sidechain

     ],
     [ // 221 : (' ', 481, ' ') C backbone
      [ 144.19919, 884, 885, 0, 0, 1,     // 481_C_N:481_C_CA:481_C_C:481_C_O [ 481_C_N:481_C_CA:481_C_C -- 481_C_CA:481_C_C:481_C_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -32.81065, 884, 886, 0, 0, 1,     // 481_C_N:481_C_CA:481_C_C:482_K_N [ 481_C_N:481_C_CA:481_C_C -- 481_C_CA:481_C_C:482_K_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 171.34516, 886, 887, 0, 0, 1,     // 481_C_CA:481_C_C:482_K_N:482_K_CA [ 481_C_CA:481_C_C:482_K_N -- 481_C_C:482_K_N:482_K_CA ] 
        [ [ 0.39524049096030783, 0.5418643993628363, 0.7417330564351352, 0.0 ], [ -0.2548190766840966, 0.8404659259619944, -0.4782094368106006, 0.0 ], [ -0.8825260293403925, 0.0, 0.4702635511462489, 15.380692577667626 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -74.47099, 887, 888, 0, 0, 1,     // 481_C_C:482_K_N:482_K_CA:482_K_C [ 481_C_C:482_K_N:482_K_CA -- 482_K_N:482_K_CA:482_K_C ] 
        [ [ -0.7866615062295741, -0.595170520453192, 0.1641210717735361, 9.927837661476556 ], [ 0.608366286122177, -0.7925497830863082, 0.04189633921663227, -6.400666136762575 ], [ 0.10513865381736195, 0.13280396422291854, 0.9855500852621092, 21.675005812576977 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 481, ' ') C sidechain

     ],
     [ // 222 : (' ', 482, ' ') K backbone
      [ 174.65498, 888, 889, 0, 0, 1,     // 482_K_N:482_K_CA:482_K_C:482_K_O [ 482_K_N:482_K_CA:482_K_C -- 482_K_CA:482_K_C:482_K_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  -5.44441, 888, 890, 0, 0, 1,     // 482_K_N:482_K_CA:482_K_C:483_N_N [ 482_K_N:482_K_CA:482_K_C -- 482_K_CA:482_K_C:483_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -174.15589, 890, 891, 0, 0, 1,     // 482_K_CA:482_K_C:483_N_N:483_N_CA [ 482_K_CA:482_K_C:483_N_N -- 482_K_C:483_N_N:483_N_CA ] 
        [ [ 0.49172377023843794, 0.09487998711594653, 0.8655665900595827, 0.0 ], [ -0.04686617149175192, 0.9954887181906572, -0.08249711464547035, 0.0 ], [ -0.8694891004217371, 0.0, 0.49395212748585104, 15.462397829898158 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  54.88712, 891, 892, 0, 0, 1,     // 482_K_C:483_N_N:483_N_CA:483_N_C [ 482_K_C:483_N_N:483_N_CA -- 483_N_N:483_N_CA:483_N_C ] 
        [ [ -0.993170804949696, -0.04431844869769944, 0.10792417384720818, 11.59428947785206 ], [ 0.034312038634610256, -0.9950867935566327, -0.09287064818292526, -1.1050512338066214 ], [ 0.11150980315788403, -0.08853331798938209, 0.9898118080753942, 22.078900846584155 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 482, ' ') K sidechain

     ],
     [ // 223 : (' ', 483, ' ') N backbone
      [ -150.31399, 892, 893, 0, 0, 1,     // 483_N_N:483_N_CA:483_N_C:483_N_O [ 483_N_N:483_N_CA:483_N_C -- 483_N_CA:483_N_C:483_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [  30.32662, 892, 894, 0, 0, 1,     // 483_N_N:483_N_CA:483_N_C:484_V_N [ 483_N_N:483_N_CA:483_N_C -- 483_N_CA:483_N_C:484_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -175.24697, 894, 895, 0, 0, 1,     // 483_N_CA:483_N_C:484_V_N:484_V_CA [ 483_N_CA:483_N_C:484_V_N -- 483_N_C:484_V_N:484_V_CA ] 
        [ [ 0.41284113900887454, -0.5049287098125657, 0.758029809406512, 0.0 ], [ 0.24150225929389257, 0.8631610498667197, 0.44342943153215214, 0.0 ], [ -0.8782020568740434, 0.0, 0.4782898151771574, 15.463730682865414 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.45800, 895, 896, 0, 0, 1,     // 483_N_C:484_V_N:484_V_CA:484_V_C [ 483_N_C:484_V_N:484_V_CA -- 484_V_N:484_V_CA:484_V_C ] 
        [ [ -0.8375929044134534, 0.5374007236323011, 0.09817631443334454, 10.162493218997135 ], [ -0.5418012014565179, -0.8401816390644344, -0.02337245128937304, 5.944817123455116 ], [ 0.06992556454194292, -0.07276864447328688, 0.9948945370264182, 21.87590200269623 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 483, ' ') N sidechain

     ],
     [ // 224 : (' ', 484, ' ') V backbone
      [ 147.39526, 896, 897, 0, 0, 1,     // 484_V_N:484_V_CA:484_V_C:484_V_O [ 484_V_N:484_V_CA:484_V_C -- 484_V_CA:484_V_C:484_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.03447, 896, 898, 0, 0, 1,     // 484_V_N:484_V_CA:484_V_C:485_V_N [ 484_V_N:484_V_CA:484_V_C -- 484_V_CA:484_V_C:485_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -174.04680, 898, 899, 0, 0, 1,     // 484_V_CA:484_V_C:485_V_N:485_V_CA [ 484_V_CA:484_V_C:485_V_N -- 484_V_C:485_V_N:485_V_CA ] 
        [ [ 0.3868946420329467, 0.5451434815316153, 0.7437278537945097, 0.0 ], [ -0.25158336424212746, 0.8383427607748453, -0.4836188853767227, 0.0 ], [ -0.8871405451239455, 0.0, 0.4614993534092859, 15.322082228282483 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -145.70374, 899, 900, 0, 0, 1,     // 484_V_C:485_V_N:485_V_CA:485_V_C [ 484_V_C:485_V_N:485_V_CA -- 485_V_N:485_V_CA:485_V_C ] 
        [ [ -0.8628544542141521, -0.50207629003726, 0.05832314999396812, 9.949754768596328 ], [ 0.49091500873387883, -0.8599147749523393, -0.13981786015557307, -6.46995979296723 ], [ 0.12035217090943318, -0.09201075372525146, 0.9884580801209011, 21.496122058659566 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 484, ' ') V sidechain

     ],
     [ // 225 : (' ', 485, ' ') V backbone
      [ -42.58571, 900, 901, 0, 0, 1,     // 485_V_N:485_V_CA:485_V_C:485_V_O [ 485_V_N:485_V_CA:485_V_C -- 485_V_CA:485_V_C:485_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 139.79140, 900, 902, 0, 0, 1,     // 485_V_N:485_V_CA:485_V_C:486_P_N [ 485_V_N:485_V_CA:485_V_C -- 485_V_CA:485_V_C:486_P_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 170.93963, 902, 903, 0, 0, 1,     // 485_V_CA:485_V_C:486_P_N:486_P_CA [ 485_V_CA:485_V_C:486_P_N -- 485_V_C:486_P_N:486_P_CA ] 
        [ [ -0.37134713120797164, -0.6455722756940178, -0.6673363057700791, 0.0 ], [ 0.31390817089210243, -0.7636991795564841, 0.564114548112285, 0.0 ], [ -0.8738209017818148, 0.0, 0.4862479116759433, 15.32234159193934 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.45959, 903, 904, 0, 0, 1,     // 485_V_C:486_P_N:486_P_CA:486_P_C [ 485_V_C:486_P_N:486_P_CA -- 486_P_N:486_P_CA:486_P_C ] 
        [ [ 0.7053411545220918, 0.6959953800633444, -0.13447782965204128, -9.009454693654439 ], [ -0.7068994304402926, 0.7047377360245569, -0.06031516116344153, 7.615896841374657 ], [ 0.05279252769694422, 0.13760506659808253, 0.9890792661186992, 21.886990468749957 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 485, ' ') V sidechain

     ],
     [ // 226 : (' ', 486, ' ') P backbone
      [ -52.23875, 904, 905, 0, 0, 1,     // 486_P_N:486_P_CA:486_P_C:486_P_O [ 486_P_N:486_P_CA:486_P_C -- 486_P_CA:486_P_C:486_P_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 129.36942, 904, 906, 0, 0, 1,     // 486_P_N:486_P_CA:486_P_C:487_V_N [ 486_P_N:486_P_CA:486_P_C -- 486_P_CA:486_P_C:487_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.52226, 906, 907, 0, 0, 1,     // 486_P_CA:486_P_C:487_V_N:487_V_CA [ 486_P_CA:486_P_C:487_V_N -- 486_P_C:487_V_N:487_V_CA ] 
        [ [ -0.2768496642351515, -0.7730721931000929, -0.5707132797371167, 0.0 ], [ 0.3374092511367666, -0.6343180466102258, 0.6955541768920043, 0.0 ], [ -0.899727325727195, 0.0, 0.43645244797089844, 15.34227478672637 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -97.80446, 907, 908, 0, 0, 1,     // 486_P_C:487_V_N:487_V_CA:487_V_C [ 486_P_C:487_V_N:487_V_CA -- 487_V_N:487_V_CA:487_V_C ] 
        [ [ 0.5707036638766099, 0.8025270996876303, -0.1739183207852864, -7.591958757031779 ], [ -0.8099535651411534, 0.5850126868452944, 0.041657874948021624, 9.252664712266121 ], [ 0.135175997695408, 0.11709146210125589, 0.9838785693112952, 21.148218159775098 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 486, ' ') P sidechain

     ],
     [ // 227 : (' ', 487, ' ') V backbone
      [ -60.23200, 908, 909, 0, 0, 1,     // 487_V_N:487_V_CA:487_V_C:487_V_O [ 487_V_N:487_V_CA:487_V_C -- 487_V_CA:487_V_C:487_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 118.24713, 908, 910, 0, 0, 1,     // 487_V_N:487_V_CA:487_V_C:488_Y_N [ 487_V_N:487_V_CA:487_V_C -- 487_V_CA:487_V_C:488_Y_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 176.44761, 910, 911, 0, 0, 1,     // 487_V_CA:487_V_C:488_Y_N:488_Y_CA [ 487_V_CA:487_V_C:488_Y_N -- 487_V_C:488_Y_N:488_Y_CA ] 
        [ [ -0.2152384617436074, -0.8809144380958738, -0.42149989008358785, 0.0 ], [ 0.4006263704541441, -0.47327555689495626, 0.7845436626106337, 0.0 ], [ -0.890601434920798, 0.0, 0.45478465686192127, 15.276107053261539 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.41897, 911, 912, 0, 0, 1,     // 487_V_C:488_Y_N:488_Y_CA:488_Y_C [ 487_V_C:488_Y_N:488_Y_CA -- 488_Y_N:488_Y_CA:488_Y_C ] 
        [ [ 0.4380671646220348, 0.8925582551897908, -0.10694353824597087, -5.622542633752022 ], [ -0.8895993161948762, 0.4475428711604148, 0.09120545542400967, 10.46531753589183 ], [ 0.12926800031569974, 0.055182783239418984, 0.9900730500969772, 21.342647478397467 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 487, ' ') V sidechain

     ],
     [ // 228 : (' ', 488, ' ') Y backbone
      [ -30.64411, 912, 913, 0, 0, 1,     // 488_Y_N:488_Y_CA:488_Y_C:488_Y_O [ 488_Y_N:488_Y_CA:488_Y_C -- 488_Y_CA:488_Y_C:488_Y_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 151.02253, 912, 914, 0, 0, 1,     // 488_Y_N:488_Y_CA:488_Y_C:489_D_N [ 488_Y_N:488_Y_CA:488_Y_C -- 488_Y_CA:488_Y_C:489_D_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -169.95278, 914, 915, 0, 0, 1,     // 488_Y_CA:488_Y_C:489_D_N:489_D_CA [ 488_Y_CA:488_Y_C:489_D_N -- 488_Y_C:489_D_N:489_D_CA ] 
        [ [ -0.4111561424143054, -0.48446564457123953, -0.7721681590075835, 0.0 ], [ 0.227696254131317, -0.8748102875653517, 0.4276229374395079, 0.0 ], [ -0.8826692712503105, 0.0, 0.4699946357039045, 15.352489350155315 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -52.89861, 915, 916, 0, 0, 1,     // 488_Y_C:489_D_N:489_D_CA:489_D_C [ 488_Y_C:489_D_N:489_D_CA -- 489_D_N:489_D_CA:489_D_C ] 
        [ [ 0.913856364221444, 0.4053058526058431, -0.024366194110162158, -10.31921376751265 ], [ -0.39488276018773083, 0.9011182929913334, 0.179034711000143, 5.714729948203721 ], [ 0.09452063943119848, -0.1539902200785368, 0.9835409807637306, 21.633472396187898 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 488, ' ') Y sidechain

     ],
     [ // 229 : (' ', 489, ' ') D backbone
      [ 128.17719, 916, 917, 0, 0, 1,     // 489_D_N:489_D_CA:489_D_C:489_D_O [ 489_D_N:489_D_CA:489_D_C -- 489_D_CA:489_D_C:489_D_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -46.51725, 916, 918, 0, 0, 1,     // 489_D_N:489_D_CA:489_D_C:490_L_N [ 489_D_N:489_D_CA:489_D_C -- 489_D_CA:489_D_C:490_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.46654, 918, 919, 0, 0, 1,     // 489_D_CA:489_D_C:490_L_N:490_L_CA [ 489_D_CA:489_D_C:490_L_N -- 489_D_C:490_L_N:490_L_CA ] 
        [ [ 0.32728450626568606, 0.7255816168936167, 0.6053231939918302, 0.0 ], [ -0.3450939646184301, 0.6881361182397309, -0.6382623585625967, 0.0 ], [ -0.879656187122196, 0.0, 0.47561012652974516, 15.39283862327874 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -64.33913, 919, 920, 0, 0, 1,     // 489_D_C:490_L_N:490_L_CA:490_L_C [ 489_D_C:490_L_N:490_L_CA -- 490_L_N:490_L_CA:490_L_C ] 
        [ [ -0.675517456910977, -0.7340801196979817, 0.06930038436201647, 8.110213121936185 ], [ 0.7333462959773667, -0.6786547258544914, -0.040385309849517624, -8.551537107831646 ], [ 0.07667708643918915, 0.023540198375584628, 0.9967780512609792, 21.765136145318394 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 489, ' ') D sidechain

     ],
     [ // 230 : (' ', 490, ' ') L backbone
      [ 130.05704, 920, 921, 0, 0, 1,     // 490_L_N:490_L_CA:490_L_C:490_L_O [ 490_L_N:490_L_CA:490_L_C -- 490_L_CA:490_L_C:490_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -48.52986, 920, 922, 0, 0, 1,     // 490_L_N:490_L_CA:490_L_C:491_L_N [ 490_L_N:490_L_CA:490_L_C -- 490_L_CA:490_L_C:491_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.82815, 922, 923, 0, 0, 1,     // 490_L_CA:490_L_C:491_L_N:491_L_CA [ 490_L_CA:490_L_C:491_L_N -- 490_L_C:491_L_N:491_L_CA ] 
        [ [ 0.3107909669144789, 0.7493009286927684, 0.5847709749504655, 0.0 ], [ -0.35165438117990966, 0.6622296567356031, -0.6616578254179573, 0.0 ], [ -0.8830335050729038, 0.0, 0.46930994973328866, 15.35935378888492 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -60.66488, 923, 924, 0, 0, 1,     // 490_L_C:491_L_N:491_L_CA:491_L_C [ 490_L_C:491_L_N:491_L_CA -- 491_L_N:491_L_CA:491_L_C ] 
        [ [ -0.6303758280740338, -0.7699250353429862, 0.09920562147417622, 7.812735980091511 ], [ 0.7723505605760381, -0.6348927239807082, -0.019642826024677932, -8.839970040561376 ], [ 0.07810843077322507, 0.06423915463685372, 0.9948730592661994, 21.629491820969687 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 490, ' ') L sidechain

     ],
     [ // 231 : (' ', 491, ' ') L backbone
      [ 142.69391, 924, 925, 0, 0, 1,     // 491_L_N:491_L_CA:491_L_C:491_L_O [ 491_L_N:491_L_CA:491_L_C -- 491_L_CA:491_L_C:491_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.30190, 924, 926, 0, 0, 1,     // 491_L_N:491_L_CA:491_L_C:492_L_N [ 491_L_N:491_L_CA:491_L_C -- 491_L_CA:491_L_C:492_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.69043, 926, 927, 0, 0, 1,     // 491_L_CA:491_L_C:492_L_N:492_L_CA [ 491_L_CA:491_L_C:492_L_N -- 491_L_C:492_L_N:492_L_CA ] 
        [ [ 0.38318256481077323, 0.5778846173865764, 0.720569560148781, 0.0 ], [ -0.2713274065156328, 0.8161184772984681, -0.5102274683761453, 0.0 ], [ -0.8829227375589204, 0.0, 0.46951830582146836, 15.414831032187875 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -69.24617, 927, 928, 0, 0, 1,     // 491_L_C:492_L_N:492_L_CA:492_L_C [ 491_L_C:492_L_N:492_L_CA -- 492_L_N:492_L_CA:492_L_C ] 
        [ [ -0.8015412543195892, -0.5864910244815533, 0.11644696572413456, 9.614881213308308 ], [ 0.5837996883758315, -0.8097043319790185, -0.05963906963303186, -6.80819281235073 ], [ 0.12926539164084758, 0.020178527621947218, 0.9914047032100204, 21.67982332108649 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 491, ' ') L sidechain

     ],
     [ // 232 : (' ', 492, ' ') L backbone
      [ 133.79916, 928, 929, 0, 0, 1,     // 492_L_N:492_L_CA:492_L_C:492_L_O [ 492_L_N:492_L_CA:492_L_C -- 492_L_CA:492_L_C:492_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -43.85207, 928, 930, 0, 0, 1,     // 492_L_N:492_L_CA:492_L_C:493_E_N [ 492_L_N:492_L_CA:492_L_C -- 492_L_CA:492_L_C:493_E_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 172.72318, 930, 931, 0, 0, 1,     // 492_L_CA:492_L_C:493_E_N:493_E_CA [ 492_L_CA:492_L_C:493_E_N -- 492_L_C:493_E_N:493_E_CA ] 
        [ [ 0.3456580225794338, 0.6927988138653742, 0.6328904604536288, 0.0 ], [ -0.33207765950323254, 0.721130919810495, -0.6080251841429963, 0.0 ], [ -0.8776360062607567, 0.0, 0.479327696377613, 15.364629080628394 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -58.23674, 931, 932, 0, 0, 1,     // 492_L_C:493_E_N:493_E_CA:493_E_C [ 492_L_C:493_E_N:493_E_CA -- 493_E_N:493_E_CA:493_E_C ] 
        [ [ -0.6663866842076531, -0.7310010579232363, 0.14685448725128833, 8.44303546118666 ], [ 0.7393950621508859, -0.6732606795324908, 0.0038729139335007996, -8.111321803356521 ], [ 0.09604024769656316, 0.11116434100267003, 0.9891505244965626, 21.75907012678415 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 492, ' ') L sidechain

     ],
     [ // 233 : (' ', 493, ' ') E backbone
      [ 122.56102, 932, 933, 0, 0, 1,     // 493_E_N:493_E_CA:493_E_C:493_E_O [ 493_E_N:493_E_CA:493_E_C -- 493_E_CA:493_E_C:493_E_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -52.16597, 932, 934, 0, 0, 1,     // 493_E_N:493_E_CA:493_E_C:494_M_N [ 493_E_N:493_E_CA:493_E_C -- 493_E_CA:493_E_C:494_M_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -179.48793, 934, 935, 0, 0, 1,     // 493_E_CA:493_E_C:494_M_N:494_M_CA [ 493_E_CA:493_E_C:494_M_N -- 493_E_C:494_M_N:494_M_CA ] 
        [ [ 0.27881182268329063, 0.7897908101024884, 0.5463464503496627, 0.0 ], [ -0.35900151012238246, 0.6133762925624491, -0.7034823661273919, 0.0 ], [ -0.8907198680066984, 0.0, 0.45455265562762887, 15.363828560875614 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -61.84716, 935, 936, 0, 0, 1,     // 493_E_C:494_M_N:494_M_CA:494_M_C [ 493_E_C:494_M_N:494_M_CA -- 494_M_N:494_M_CA:494_M_C ] 
        [ [ -0.6123554122877093, -0.7872674662636001, 0.07232417026721315, 7.304466652494215 ], [ 0.7802538180092575, -0.6165602696553083, -0.10515423607515706, -9.40532052639231 ], [ 0.12737671892432273, -0.007960555600257234, 0.9918224644716459, 21.44104332581044 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 493, ' ') E sidechain

     ],
     [ // 234 : (' ', 494, ' ') M backbone
      [ 142.48823, 936, 937, 0, 0, 1,     // 494_M_N:494_M_CA:494_M_C:494_M_O [ 494_M_N:494_M_CA:494_M_C -- 494_M_CA:494_M_C:494_M_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -35.05048, 936, 938, 0, 0, 1,     // 494_M_N:494_M_CA:494_M_C:495_L_N [ 494_M_N:494_M_CA:494_M_C -- 494_M_CA:494_M_C:495_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.59765, 938, 939, 0, 0, 1,     // 494_M_CA:494_M_C:495_L_N:495_L_CA [ 494_M_CA:494_M_C:495_L_N -- 494_M_C:495_L_N:495_L_CA ] 
        [ [ 0.3811954082731362, 0.5742978711626356, 0.7244805144995563, 0.0 ], [ -0.26741668572332156, 0.818646416457108, -0.5082384883287466, 0.0 ], [ -0.8849736588781789, 0.0, 0.4656410882769783, 15.395807833908698 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.02919, 939, 940, 0, 0, 1,     // 494_M_C:495_L_N:495_L_CA:495_L_C [ 494_M_C:495_L_N:495_L_CA -- 495_L_N:495_L_CA:495_L_C ] 
        [ [ -0.7862406027558282, -0.6018640185165101, 0.13994791100016463, 9.67379135430171 ], [ 0.605309268289255, -0.7957041792425096, -0.021343590586655047, -6.786370365963806 ], [ 0.12420307685914515, 0.06793057007829416, 0.9899288021608217, 21.613386697424467 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 494, ' ') M sidechain

     ],
     [ // 235 : (' ', 495, ' ') L backbone
      [ 135.44555, 940, 941, 0, 0, 1,     // 495_L_N:495_L_CA:495_L_C:495_L_O [ 495_L_N:495_L_CA:495_L_C -- 495_L_CA:495_L_C:495_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -42.06058, 940, 942, 0, 0, 1,     // 495_L_N:495_L_CA:495_L_C:496_N_N [ 495_L_N:495_L_CA:495_L_C -- 495_L_CA:495_L_C:496_N_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.31682, 942, 943, 0, 0, 1,     // 495_L_CA:495_L_C:496_N_N:496_N_CA [ 495_L_CA:495_L_C:496_N_N -- 495_L_C:496_N_N:496_N_CA ] 
        [ [ 0.34966803824849174, 0.6699159259351991, 0.6549388637161765, 0.0 ], [ -0.31551255871603767, 0.7424369684884939, -0.5909646124042235, 0.0 ], [ -0.8821474300364482, 0.0, 0.47097336621096714, 15.430287792145199 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -65.71623, 943, 944, 0, 0, 1,     // 495_L_C:496_N_N:496_N_CA:496_N_C [ 495_L_C:496_N_N:496_N_CA -- 496_N_N:496_N_CA:496_N_C ] 
        [ [ -0.7276066437570701, -0.6798976144505713, 0.09125681248539512, 8.757062948651607 ], [ 0.6785883592612001, -0.7328491285780616, -0.049497408190411964, -7.901675405067984 ], [ 0.10053064525688064, 0.025911167607771908, 0.994596501480594, 21.727583071960975 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 495, ' ') L sidechain

     ],
     [ // 236 : (' ', 496, ' ') N backbone
      [ 143.76381, 944, 945, 0, 0, 1,     // 496_N_N:496_N_CA:496_N_C:496_N_O [ 496_N_N:496_N_CA:496_N_C -- 496_N_CA:496_N_C:496_N_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.67188, 944, 946, 0, 0, 1,     // 496_N_N:496_N_CA:496_N_C:497_A_N [ 496_N_N:496_N_CA:496_N_C -- 496_N_CA:496_N_C:497_A_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.69752, 946, 947, 0, 0, 1,     // 496_N_CA:496_N_C:497_A_N:497_A_CA [ 496_N_CA:496_N_C:497_A_N -- 496_N_C:497_A_N:497_A_CA ] 
        [ [ 0.39099445286528, 0.5544360101147121, 0.7346591376391226, 0.0 ], [ -0.260483704199852, 0.8322263578426718, -0.4894359295742857, 0.0 ], [ -0.8827636023731973, 0.0, 0.469817435101227, 15.40750501581741 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -57.23328, 947, 948, 0, 0, 1,     // 496_N_C:497_A_N:497_A_CA:497_A_C [ 496_N_C:497_A_N:497_A_CA -- 497_A_N:497_A_CA:497_A_C ] 
        [ [ -0.8151406885728416, -0.5696965951535597, 0.10484010350736417, 9.827718802671116 ], [ 0.5700206906889537, -0.8210895328465528, -0.029805892644204638, -6.547306691423471 ], [ 0.1030634271673618, 0.035465032359643954, 0.9940423338370692, 21.692370205958298 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 496, ' ') N sidechain

     ],
     [ // 237 : (' ', 497, ' ') A backbone
      [ 142.81343, 948, 949, 0, 0, 1,     // 497_A_N:497_A_CA:497_A_C:497_A_O [ 497_A_N:497_A_CA:497_A_C -- 497_A_CA:497_A_C:497_A_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -33.84778, 948, 950, 0, 0, 1,     // 497_A_N:497_A_CA:497_A_C:498_H_N [ 497_A_N:497_A_CA:497_A_C -- 497_A_CA:497_A_C:498_H_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.89804, 950, 951, 0, 0, 1,     // 497_A_CA:497_A_C:498_H_N:498_H_CA [ 497_A_CA:497_A_C:498_H_N -- 497_A_C:498_H_N:498_H_CA ] 
        [ [ 0.3979924849520451, 0.556988385133947, 0.7289485034949822, 0.0 ], [ -0.2669136399090067, 0.8305202820075365, -0.4888692770103754, 0.0 ], [ -0.8777010258352336, 0.0, 0.4792086281023942, 15.435501080428491 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -72.24665, 951, 952, 0, 0, 1,     // 497_A_C:498_H_N:498_H_CA:498_H_C [ 497_A_C:498_H_N:498_H_CA -- 498_H_N:498_H_CA:498_H_C ] 
        [ [ -0.8004760938751543, -0.5901744508032448, 0.1045568780782195, 9.763278883031935 ], [ 0.5949843444531362, -0.8034935685195039, 0.019791796367810872, -6.547742489235172 ], [ 0.07233016652854485, 0.07805256540871364, 0.9943219519064594, 21.853852249114325 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 497, ' ') A sidechain

     ],
     [ // 238 : (' ', 498, ' ') H backbone
      [ 135.51401, 952, 953, 0, 0, 1,     // 498_H_N:498_H_CA:498_H_C:498_H_O [ 498_H_N:498_H_CA:498_H_C -- 498_H_CA:498_H_C:498_H_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -40.88201, 952, 954, 0, 0, 1,     // 498_H_N:498_H_CA:498_H_C:499_V_N [ 498_H_N:498_H_CA:498_H_C -- 498_H_CA:498_H_C:499_V_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 178.69071, 954, 955, 0, 0, 1,     // 498_H_CA:498_H_C:499_V_N:499_V_CA [ 498_H_CA:498_H_C:499_V_N -- 498_H_C:499_V_N:499_V_CA ] 
        [ [ 0.3668155137430027, 0.6545034268318157, 0.661113941119734, 0.0 ], [ -0.31754400076603756, 0.7560590349075991, -0.5723115788729848, 0.0 ], [ -0.8744210578748408, 0.0, 0.485167819980926, 15.365263574514255 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -59.70999, 955, 956, 0, 0, 1,     // 498_H_C:499_V_N:499_V_CA:499_V_C [ 498_H_C:499_V_N:499_V_CA -- 499_V_N:499_V_CA:499_V_C ] 
        [ [ -0.7453175427822755, -0.6627140556998548, 0.07288237646263154, 8.855615490630738 ], [ 0.6619298541281124, -0.7486059582876733, -0.037920804713640684, -7.666108620778647 ], [ 0.07969083156127024, 0.019979979830922676, 0.996619371561194, 21.86408230294961 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 498, ' ') H sidechain

     ],
     [ // 239 : (' ', 499, ' ') V backbone
      [ 130.00449, 956, 957, 0, 0, 1,     // 499_V_N:499_V_CA:499_V_C:499_V_O [ 499_V_N:499_V_CA:499_V_C -- 499_V_CA:499_V_C:499_V_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -45.02060, 956, 958, 0, 0, 1,     // 499_V_N:499_V_CA:499_V_C:500_L_N [ 499_V_N:499_V_CA:499_V_C -- 499_V_CA:499_V_C:500_L_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 177.32300, 958, 959, 0, 0, 1,     // 499_V_CA:499_V_C:500_L_N:500_L_CA [ 499_V_CA:499_V_C:500_L_N -- 499_V_C:500_L_N:500_L_CA ] 
        [ [ 0.31698618576154297, 0.7073609991353691, 0.63179124316388, 0.0 ], [ -0.31721423354242106, 0.7068524718088014, -0.6322457695674064, 0.0 ], [ -0.893809201157855, 0.0, 0.4484474461133178, 15.309208946572186 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.23975, 959, 960, 0, 0, 1,     // 499_V_C:500_L_N:500_L_CA:500_L_C [ 499_V_C:500_L_N:500_L_CA -- 500_L_N:500_L_CA:500_L_C ] 
        [ [ -0.681695277832612, -0.7213940586036247, 0.12199246039042028, 8.446108996649134 ], [ 0.7194405102616611, -0.6912654333233296, -0.06750890968420502, -8.45218533845915 ], [ 0.1330296973429572, 0.0417458130080086, 0.9902324912469503, 21.304284139723524 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 499, ' ') V sidechain

     ],
     [ // 240 : (' ', 500, ' ') L backbone
      [ 155.60320, 960, 961, 0, 0, 1,     // 500_L_N:500_L_CA:500_L_C:500_L_O [ 500_L_N:500_L_CA:500_L_C -- 500_L_CA:500_L_C:500_L_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -24.17262, 960, 962, 0, 0, 1,     // 500_L_N:500_L_CA:500_L_C:501_R_N [ 500_L_N:500_L_CA:500_L_C -- 500_L_CA:500_L_C:501_R_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 174.79699, 962, 963, 0, 0, 1,     // 500_L_CA:500_L_C:501_R_N:501_R_CA [ 500_L_CA:500_L_C:501_R_N -- 500_L_C:501_R_N:501_R_CA ] 
        [ [ 0.4330708736807463, 0.4094871436071387, 0.802975651928414, 0.0 ], [ -0.19438108817903293, 0.912315887848319, -0.36041047895124867, 0.0 ], [ -0.8801509023614814, 0.0, 0.47469399519297684, 15.382511447138148 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -84.17649, 963, 964, 0, 0, 1,     // 500_L_C:501_R_N:501_R_CA:501_R_C [ 500_L_C:501_R_N:501_R_CA -- 501_R_N:501_R_CA:501_R_C ] 
        [ [ -0.885241076435617, -0.44707290652695775, 0.12835128686790886, 10.750531384574606 ], [ 0.45348524924883865, -0.8909293773674685, 0.024412563514176932, -4.825307163411212 ], [ 0.10343773636748232, 0.07981641932054569, 0.9914282494976719, 21.737888077967995 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 500, ' ') L sidechain

     ],
     [ // 241 : (' ', 501, ' ') R backbone
      [ 137.36228, 964, 965, 0, 0, 1,     // 501_R_N:501_R_CA:501_R_C:501_R_O [ 501_R_N:501_R_CA:501_R_C -- 501_R_CA:501_R_C:501_R_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -38.74354, 964, 966, 0, 0, 1,     // 501_R_N:501_R_CA:501_R_C:502_G_N [ 501_R_N:501_R_CA:501_R_C -- 501_R_CA:501_R_C:502_G_N ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ 175.31235, 966, 967, 0, 0, 1,     // 501_R_CA:501_R_C:502_G_N:502_G_CA [ 501_R_CA:501_R_C:502_G_N -- 501_R_C:502_G_N:502_G_CA ] 
        [ [ 0.3737808775917437, 0.6258355833297885, 0.6845565558703016, 0.0 ], [ -0.2999216231729614, 0.7799550132143926, -0.5492879002081084, 0.0 ], [ -0.8776872310225563, 0.0, 0.479233893313023, 15.331680392882037 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
      [ -63.15659, 967, 968, 0, 0, 1,     // 501_R_C:502_G_N:502_G_CA:502_G_C [ 501_R_C:502_G_N:502_G_CA -- 502_G_N:502_G_CA:502_G_C ] 
        [ [ -0.7474494457972626, -0.6542889211210057, 0.11500145075465634, 9.14920384552518 ], [ 0.65818746229412, -0.7528353928973437, -0.005304307676658429, -7.341317420436235 ], [ 0.09004771210970093, 0.07172781119909295, 0.9933511617974762, 21.73671521725999 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 501, ' ') R sidechain

     ],
     [ // 242 : (' ', 502, ' ') G backbone
      [ 119.43778, 968, 969, 0, 0, 1,     // 502_G_N:502_G_CA:502_G_C:502_G_O [ 502_G_N:502_G_CA:502_G_C -- 502_G_CA:502_G_C:502_G_O ] 
        [ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
       // (' ', 502, ' ') G sidechain
   ],
  ],
   [  //hedra
     [  14.68326, 113.29741,  15.27609, "Nsb", "Csb", "Cdb", 1, 1, 1, StdBond, StdBond, "L", 260, "NCAC" ], // (260_L_N, 260_L_CA, 260_L_C)
     [  15.27609, 118.98079,  12.24945, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 260, "CACO" ], // (260_L_CA, 260_L_C, 260_L_O)
     [  15.27609, 118.07139,  13.35942, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 260, "CACN" ], // (260_L_CA, 260_L_C, 261_D_N)
     [  13.35942, 123.33452,  14.73105, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 260, "CNCA" ], // (260_L_C, 261_D_N, 261_D_CA)
     [  14.73105, 112.23529,  15.45191, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 261, "NCAC" ], // (261_D_N, 261_D_CA, 261_D_C)
     [  15.45191, 120.76051,  12.30376, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 261, "CACO" ], // (261_D_CA, 261_D_C, 261_D_O)
     [  15.45191, 117.64965,  13.39317, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 261, "CACN" ], // (261_D_CA, 261_D_C, 262_A_N)
     [  13.39317, 125.16877,  14.59768, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 261, "CNCA" ], // (261_D_C, 262_A_N, 262_A_CA)
     [  14.59768, 114.78051,  15.45449, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 262, "NCAC" ], // (262_A_N, 262_A_CA, 262_A_C)
     [  15.45449, 118.77973,  12.32587, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 262, "CACO" ], // (262_A_CA, 262_A_C, 262_A_O)
     [  15.45449, 119.68584,  13.38862, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 262, "CACN" ], // (262_A_CA, 262_A_C, 263_L_N)
     [  13.38862, 123.27491,  14.65732, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 262, "CNCA" ], // (262_A_C, 263_L_N, 263_L_CA)
     [  14.65732, 110.91721,  15.27350, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 263, "NCAC" ], // (263_L_N, 263_L_CA, 263_L_C)
     [  15.27350, 121.97273,  12.29262, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 263, "CACO" ], // (263_L_CA, 263_L_C, 263_L_O)
     [  15.27350, 116.20068,  13.31295, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 263, "CACN" ], // (263_L_CA, 263_L_C, 264_S_N)
     [  13.31295, 121.36786,  14.61968, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 263, "CNCA" ], // (263_L_C, 264_S_N, 264_S_CA)
     [  14.61968, 111.28059,  15.41347, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 264, "NCAC" ], // (264_S_N, 264_S_CA, 264_S_C)
     [  15.41347, 119.91865,  12.36603, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 264, "CACO" ], // (264_S_CA, 264_S_C, 264_S_O)
     [  15.41347, 118.45193,  13.48562, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 264, "CACN" ], // (264_S_CA, 264_S_C, 265_P_N)
     [  13.48562, 123.55939,  14.65425, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 264, "CNCA" ], // (264_S_C, 265_P_N, 265_P_CA)
     [  14.65425, 113.82490,  15.38255, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 265, "NCAC" ], // (265_P_N, 265_P_CA, 265_P_C)
     [  15.38255, 120.54137,  12.27709, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 265, "CACO" ], // (265_P_CA, 265_P_C, 265_P_O)
     [  15.38255, 117.57370,  13.36235, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 265, "CACN" ], // (265_P_CA, 265_P_C, 266_E_N)
     [  13.36235, 124.24393,  14.66078, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 265, "CNCA" ], // (265_P_C, 266_E_N, 266_E_CA)
     [  14.66078, 112.62854,  15.34353, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 266, "NCAC" ], // (266_E_N, 266_E_CA, 266_E_C)
     [  15.34353, 120.09308,  12.30815, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 266, "CACO" ], // (266_E_CA, 266_E_C, 266_E_O)
     [  15.34353, 117.44581,  13.36765, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 266, "CACN" ], // (266_E_CA, 266_E_C, 267_Q_N)
     [  13.36765, 122.98840,  14.64665, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 266, "CNCA" ], // (266_E_C, 267_Q_N, 267_Q_CA)
     [  14.64665, 111.37053,  15.29799, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 267, "NCAC" ], // (267_Q_N, 267_Q_CA, 267_Q_C)
     [  15.29799, 118.91765,  12.26128, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 267, "CACO" ], // (267_Q_CA, 267_Q_C, 267_Q_O)
     [  15.29799, 118.63433,  13.39028, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 267, "CACN" ], // (267_Q_CA, 267_Q_C, 268_L_N)
     [  13.39028, 122.21058,  14.64407, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 267, "CNCA" ], // (267_Q_C, 268_L_N, 268_L_CA)
     [  14.64407, 111.64716,  15.38739, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 268, "NCAC" ], // (268_L_N, 268_L_CA, 268_L_C)
     [  15.38739, 119.60654,  12.27864, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 268, "CACO" ], // (268_L_CA, 268_L_C, 268_L_O)
     [  15.38739, 118.94911,  13.39463, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 268, "CACN" ], // (268_L_CA, 268_L_C, 269_V_N)
     [  13.39463, 122.75444,  14.77961, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 268, "CNCA" ], // (268_L_C, 269_V_N, 269_V_CA)
     [  14.77961, 112.19036,  15.32247, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 269, "NCAC" ], // (269_V_N, 269_V_CA, 269_V_C)
     [  15.32247, 119.28863,  12.30882, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 269, "CACO" ], // (269_V_CA, 269_V_C, 269_V_O)
     [  15.32247, 118.72699,  13.36901, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 269, "CACN" ], // (269_V_CA, 269_V_C, 270_L_N)
     [  13.36901, 122.71381,  14.68012, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 269, "CNCA" ], // (269_V_C, 270_L_N, 270_L_CA)
     [  14.68012, 111.48548,  15.35869, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 270, "NCAC" ], // (270_L_N, 270_L_CA, 270_L_C)
     [  15.35869, 119.26906,  12.27864, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 270, "CACO" ], // (270_L_CA, 270_L_C, 270_L_O)
     [  15.35869, 118.74757,  13.36507, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 270, "CACN" ], // (270_L_CA, 270_L_C, 271_T_N)
     [  13.36507, 122.55960,  14.71706, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 270, "CNCA" ], // (270_L_C, 271_T_N, 271_T_CA)
     [  14.71706, 111.98936,  15.31127, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 271, "NCAC" ], // (271_T_N, 271_T_CA, 271_T_C)
     [  15.31127, 119.21319,  12.26146, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 271, "CACO" ], // (271_T_CA, 271_T_C, 271_T_O)
     [  15.31127, 119.45134,  13.40510, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 271, "CACN" ], // (271_T_CA, 271_T_C, 272_L_N)
     [  13.40510, 122.42463,  14.70482, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 271, "CNCA" ], // (271_T_C, 272_L_N, 272_L_CA)
     [  14.70482, 112.55606,  15.40123, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 272, "NCAC" ], // (272_L_N, 272_L_CA, 272_L_C)
     [  15.40123, 119.24477,  12.28414, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 272, "CACO" ], // (272_L_CA, 272_L_C, 272_L_O)
     [  15.40123, 118.48976,  13.37204, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 272, "CACN" ], // (272_L_CA, 272_L_C, 273_L_N)
     [  13.37204, 123.79683,  14.66795, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 272, "CNCA" ], // (272_L_C, 273_L_N, 273_L_CA)
     [  14.66795, 111.67707,  15.38741, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 273, "NCAC" ], // (273_L_N, 273_L_CA, 273_L_C)
     [  15.38741, 119.82771,  12.27664, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 273, "CACO" ], // (273_L_CA, 273_L_C, 273_L_O)
     [  15.38741, 118.11320,  13.35257, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 273, "CACN" ], // (273_L_CA, 273_L_C, 274_E_N)
     [  13.35257, 124.32267,  14.76648, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 273, "CNCA" ], // (273_L_C, 274_E_N, 274_E_CA)
     [  14.76648, 112.74743,  15.40882, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 274, "NCAC" ], // (274_E_N, 274_E_CA, 274_E_C)
     [  15.40882, 121.89580,  12.38311, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 274, "CACO" ], // (274_E_CA, 274_E_C, 274_E_O)
     [  15.40882, 117.09912,  13.34976, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 274, "CACN" ], // (274_E_CA, 274_E_C, 275_A_N)
     [  13.34976, 122.77154,  14.48774, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 274, "CNCA" ], // (274_E_C, 275_A_N, 275_A_CA)
     [  14.48774, 113.31481,  15.32674, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 275, "NCAC" ], // (275_A_N, 275_A_CA, 275_A_C)
     [  15.32674, 117.35669,  12.38988, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 275, "CACO" ], // (275_A_CA, 275_A_C, 275_A_O)
     [  15.32674, 118.99100,  13.41702, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 275, "CACN" ], // (275_A_CA, 275_A_C, 276_E_N)
     [  13.41702, 123.36522,  14.67400, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 275, "CNCA" ], // (275_A_C, 276_E_N, 276_E_CA)
     [  14.67400, 113.49360,  15.39805, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 276, "NCAC" ], // (276_E_N, 276_E_CA, 276_E_C)
     [  15.39805, 121.05296,  12.37984, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 276, "CACO" ], // (276_E_CA, 276_E_C, 276_E_O)
     [  15.39805, 117.79422,  13.45747, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 276, "CACN" ], // (276_E_CA, 276_E_C, 277_P_N)
     [  13.45747, 121.35566,  14.49864, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 276, "CNCA" ], // (276_E_C, 277_P_N, 277_P_CA)
     [  14.49864, 111.50219,  15.40363, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 277, "NCAC" ], // (277_P_N, 277_P_CA, 277_P_C)
     [  15.40363, 120.15887,  12.32258, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 277, "CACO" ], // (277_P_CA, 277_P_C, 277_P_O)
     [  15.40363, 118.12937,  13.45715, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 277, "CACN" ], // (277_P_CA, 277_P_C, 278_P_N)
     [  13.45715, 122.74345,  14.59159, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 277, "CNCA" ], // (277_P_C, 278_P_N, 278_P_CA)
     [  14.59159, 111.56103,  15.34915, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 278, "NCAC" ], // (278_P_N, 278_P_CA, 278_P_C)
     [  15.34915, 121.40764,  12.31237, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 278, "CACO" ], // (278_P_CA, 278_P_C, 278_P_O)
     [  15.34915, 116.63311,  13.33839, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 278, "CACN" ], // (278_P_CA, 278_P_C, 279_H_N)
     [  13.33839, 121.92490,  14.61738, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 278, "CNCA" ], // (278_P_C, 279_H_N, 279_H_CA)
     [  14.61738, 111.23173,  15.36668, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 279, "NCAC" ], // (279_H_N, 279_H_CA, 279_H_C)
     [  15.36668, 121.05315,  12.31268, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 279, "CACO" ], // (279_H_CA, 279_H_C, 279_H_O)
     [  15.36668, 116.21011,  13.32356, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 279, "CACN" ], // (279_H_CA, 279_H_C, 280_V_N)
     [  13.32356, 125.01256,  14.67541, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 279, "CNCA" ], // (279_H_C, 280_V_N, 280_V_CA)
     [  14.67541, 110.01987,  15.29691, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 280, "NCAC" ], // (280_V_N, 280_V_CA, 280_V_C)
     [  15.29691, 119.57159,  12.25277, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 280, "CACO" ], // (280_V_CA, 280_V_C, 280_V_O)
     [  15.29691, 117.39312,  13.37795, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 280, "CACN" ], // (280_V_CA, 280_V_C, 281_L_N)
     [  13.37795, 125.51983,  14.73216, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 280, "CNCA" ], // (280_V_C, 281_L_N, 281_L_CA)
     [  14.73216, 110.58947,  15.40304, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 281, "NCAC" ], // (281_L_N, 281_L_CA, 281_L_C)
     [  15.40304, 120.34548,  12.40583, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 281, "CACO" ], // (281_L_CA, 281_L_C, 281_L_O)
     [  15.40304, 118.29557,  13.35829, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 281, "CACN" ], // (281_L_CA, 281_L_C, 282_I_N)
     [  13.35829, 123.82561,  14.84904, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 281, "CNCA" ], // (281_L_C, 282_I_N, 282_I_CA)
     [  14.84904, 112.70692,  15.38598, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 282, "NCAC" ], // (282_I_N, 282_I_CA, 282_I_C)
     [  15.38598, 119.18689,  12.41043, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 282, "CACO" ], // (282_I_CA, 282_I_C, 282_I_O)
     [  15.38598, 115.83157,  13.37292, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 282, "CACN" ], // (282_I_CA, 282_I_C, 283_S_N)
     [  13.37292, 129.02178,  14.65895, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 282, "CNCA" ], // (282_I_C, 283_S_N, 283_S_CA)
     [  14.65895, 112.12718,  15.34239, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 283, "NCAC" ], // (283_S_N, 283_S_CA, 283_S_C)
     [  15.34239, 120.12589,  12.31442, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 283, "CACO" ], // (283_S_CA, 283_S_C, 283_S_O)
     [  15.34239, 117.46496,  13.43357, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 283, "CACN" ], // (283_S_CA, 283_S_C, 284_R_N)
     [  13.43357, 124.30508,  14.78141, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 283, "CNCA" ], // (283_S_C, 284_R_N, 284_R_CA)
     [  14.78141, 110.81050,  15.39810, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 284, "NCAC" ], // (284_R_N, 284_R_CA, 284_R_C)
     [  15.39810, 118.39721,  12.32917, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 284, "CACO" ], // (284_R_CA, 284_R_C, 284_R_O)
     [  15.39810, 119.72884,  13.47277, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 284, "CACN" ], // (284_R_CA, 284_R_C, 285_P_N)
     [  13.47277, 121.86206,  14.60941, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 284, "CNCA" ], // (284_R_C, 285_P_N, 285_P_CA)
     [  14.60941, 112.45714,  15.32093, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 285, "NCAC" ], // (285_P_N, 285_P_CA, 285_P_C)
     [  15.32093, 121.17373,  12.32779, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 285, "CACO" ], // (285_P_CA, 285_P_C, 285_P_O)
     [  15.32093, 116.45197,  13.31195, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 285, "CACN" ], // (285_P_CA, 285_P_C, 286_S_N)
     [  13.31195, 122.84213,  14.62565, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 285, "CNCA" ], // (285_P_C, 286_S_N, 286_S_CA)
     [  14.62565, 110.44717,  15.29953, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 286, "NCAC" ], // (286_S_N, 286_S_CA, 286_S_C)
     [  15.29953, 119.77274,  12.30432, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 286, "CACO" ], // (286_S_CA, 286_S_C, 286_S_O)
     [  15.29953, 117.33287,  13.32040, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 286, "CACN" ], // (286_S_CA, 286_S_C, 287_A_N)
     [  13.32040, 122.64572,  14.59957, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 286, "CNCA" ], // (286_S_C, 287_A_N, 287_A_CA)
     [  14.59957, 111.11723,  15.49857, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 287, "NCAC" ], // (287_A_N, 287_A_CA, 287_A_C)
     [  15.49857, 118.84603,  12.35189, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 287, "CACO" ], // (287_A_CA, 287_A_C, 287_A_O)
     [  15.49857, 119.54505,  13.48623, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 287, "CACN" ], // (287_A_CA, 287_A_C, 288_P_N)
     [  13.48623, 127.21302,  14.66025, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 287, "CNCA" ], // (287_A_C, 288_P_N, 288_P_CA)
     [  14.66025, 112.50392,  15.36062, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 288, "NCAC" ], // (288_P_N, 288_P_CA, 288_P_C)
     [  15.36062, 121.14861,  12.28848, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 288, "CACO" ], // (288_P_CA, 288_P_C, 288_P_O)
     [  15.36062, 116.75522,  13.34419, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 288, "CACN" ], // (288_P_CA, 288_P_C, 289_F_N)
     [  13.34419, 121.98944,  14.61830, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 288, "CNCA" ], // (288_P_C, 289_F_N, 289_F_CA)
     [  14.61830, 111.77774,  15.33722, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 289, "NCAC" ], // (289_F_N, 289_F_CA, 289_F_C)
     [  15.33722, 121.36065,  12.30396, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 289, "CACO" ], // (289_F_CA, 289_F_C, 289_F_O)
     [  15.33722, 116.00817,  13.30268, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 289, "CACN" ], // (289_F_CA, 289_F_C, 290_T_N)
     [  13.30268, 124.63135,  14.70525, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 289, "CNCA" ], // (289_F_C, 290_T_N, 290_T_CA)
     [  14.70525, 110.37396,  15.29604, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 290, "NCAC" ], // (290_T_N, 290_T_CA, 290_T_C)
     [  15.29604, 119.25234,  12.35246, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 290, "CACO" ], // (290_T_CA, 290_T_C, 290_T_O)
     [  15.29604, 117.74095,  13.38591, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 290, "CACN" ], // (290_T_CA, 290_T_C, 291_E_N)
     [  13.38591, 124.49078,  14.74360, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 290, "CNCA" ], // (290_T_C, 291_E_N, 291_E_CA)
     [  14.74360, 112.36573,  15.38077, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 291, "NCAC" ], // (291_E_N, 291_E_CA, 291_E_C)
     [  15.38077, 121.11142,  12.32050, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 291, "CACO" ], // (291_E_CA, 291_E_C, 291_E_O)
     [  15.38077, 116.66969,  13.37596, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 291, "CACN" ], // (291_E_CA, 291_E_C, 292_A_N)
     [  13.37596, 123.93685,  14.59589, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 291, "CNCA" ], // (291_E_C, 292_A_N, 292_A_CA)
     [  14.59589, 112.87828,  15.40195, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 292, "NCAC" ], // (292_A_N, 292_A_CA, 292_A_C)
     [  15.40195, 119.51728,  12.29182, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 292, "CACO" ], // (292_A_CA, 292_A_C, 292_A_O)
     [  15.40195, 119.47551,  13.41064, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 292, "CACN" ], // (292_A_CA, 292_A_C, 293_S_N)
     [  13.41064, 122.24078,  14.76097, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 292, "CNCA" ], // (292_A_C, 293_S_N, 293_S_CA)
     [  14.76097, 112.77497,  15.33167, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 293, "NCAC" ], // (293_S_N, 293_S_CA, 293_S_C)
     [  15.33167, 120.35888,  12.35709, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 293, "CACO" ], // (293_S_CA, 293_S_C, 293_S_O)
     [  15.33167, 116.71881,  13.36839, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 293, "CACN" ], // (293_S_CA, 293_S_C, 294_M_N)
     [  13.36839, 124.54793,  14.56230, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 293, "CNCA" ], // (293_S_C, 294_M_N, 294_M_CA)
     [  14.56230, 109.68404,  15.25474, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 294, "NCAC" ], // (294_M_N, 294_M_CA, 294_M_C)
     [  15.25474, 120.52980,  12.30829, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 294, "CACO" ], // (294_M_CA, 294_M_C, 294_M_O)
     [  15.25474, 117.08215,  13.34440, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 294, "CACN" ], // (294_M_CA, 294_M_C, 295_M_N)
     [  13.34440, 122.75482,  14.58645, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 294, "CNCA" ], // (294_M_C, 295_M_N, 295_M_CA)
     [  14.58645, 112.24305,  15.33574, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 295, "NCAC" ], // (295_M_N, 295_M_CA, 295_M_C)
     [  15.33574, 118.67842,  12.24872, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 295, "CACO" ], // (295_M_CA, 295_M_C, 295_M_O)
     [  15.33574, 119.66522,  13.38898, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 295, "CACN" ], // (295_M_CA, 295_M_C, 296_M_N)
     [  13.38898, 123.59022,  14.71618, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 295, "CNCA" ], // (295_M_C, 296_M_N, 296_M_CA)
     [  14.71618, 112.50197,  15.33003, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 296, "NCAC" ], // (296_M_N, 296_M_CA, 296_M_C)
     [  15.33003, 120.28993,  12.29648, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 296, "CACO" ], // (296_M_CA, 296_M_C, 296_M_O)
     [  15.33003, 118.19386,  13.33229, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 296, "CACN" ], // (296_M_CA, 296_M_C, 297_S_N)
     [  13.33229, 122.41134,  14.67958, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 296, "CNCA" ], // (296_M_C, 297_S_N, 297_S_CA)
     [  14.67958, 111.46612,  15.26807, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 297, "NCAC" ], // (297_S_N, 297_S_CA, 297_S_C)
     [  15.26807, 119.88861,  12.31177, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 297, "CACO" ], // (297_S_CA, 297_S_C, 297_S_O)
     [  15.26807, 117.24236,  13.35662, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 297, "CACN" ], // (297_S_CA, 297_S_C, 298_L_N)
     [  13.35662, 124.04143,  14.62494, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 297, "CNCA" ], // (297_S_C, 298_L_N, 298_L_CA)
     [  14.62494, 112.38822,  15.38209, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 298, "NCAC" ], // (298_L_N, 298_L_CA, 298_L_C)
     [  15.38209, 118.82997,  12.30123, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 298, "CACO" ], // (298_L_CA, 298_L_C, 298_L_O)
     [  15.38209, 119.97521,  13.39319, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 298, "CACN" ], // (298_L_CA, 298_L_C, 299_T_N)
     [  13.39319, 122.36928,  14.77565, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 298, "CNCA" ], // (298_L_C, 299_T_N, 299_T_CA)
     [  14.77565, 111.86017,  15.27660, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 299, "NCAC" ], // (299_T_N, 299_T_CA, 299_T_C)
     [  15.27660, 119.56288,  12.29537, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 299, "CACO" ], // (299_T_CA, 299_T_C, 299_T_O)
     [  15.27660, 118.29311,  13.39605, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 299, "CACN" ], // (299_T_CA, 299_T_C, 300_K_N)
     [  13.39605, 124.11877,  14.74856, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 299, "CNCA" ], // (299_T_C, 300_K_N, 300_K_CA)
     [  14.74856, 112.89829,  15.30800, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 300, "NCAC" ], // (300_K_N, 300_K_CA, 300_K_C)
     [  15.30800, 119.32668,  12.29631, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 300, "CACO" ], // (300_K_CA, 300_K_C, 300_K_O)
     [  15.30800, 118.80683,  13.36448, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 300, "CACN" ], // (300_K_CA, 300_K_C, 301_L_N)
     [  13.36448, 122.34585,  14.64494, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 300, "CNCA" ], // (300_K_C, 301_L_N, 301_L_CA)
     [  14.64494, 110.68824,  15.27643, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 301, "NCAC" ], // (301_L_N, 301_L_CA, 301_L_C)
     [  15.27643, 119.04450,  12.27072, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 301, "CACO" ], // (301_L_CA, 301_L_C, 301_L_O)
     [  15.27643, 117.47058,  13.35892, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 301, "CACN" ], // (301_L_CA, 301_L_C, 302_A_N)
     [  13.35892, 123.32166,  14.58165, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 301, "CNCA" ], // (301_L_C, 302_A_N, 302_A_CA)
     [  14.58165, 111.91038,  15.35604, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 302, "NCAC" ], // (302_A_N, 302_A_CA, 302_A_C)
     [  15.35604, 119.30159,  12.24605, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 302, "CACO" ], // (302_A_CA, 302_A_C, 302_A_O)
     [  15.35604, 118.53879,  13.35517, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 302, "CACN" ], // (302_A_CA, 302_A_C, 303_D_N)
     [  13.35517, 123.80646,  14.68074, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 302, "CNCA" ], // (302_A_C, 303_D_N, 303_D_CA)
     [  14.68074, 112.56800,  15.45026, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 303, "NCAC" ], // (303_D_N, 303_D_CA, 303_D_C)
     [  15.45026, 121.15047,  12.33503, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 303, "CACO" ], // (303_D_CA, 303_D_C, 303_D_O)
     [  15.45026, 117.93192,  13.40836, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 303, "CACN" ], // (303_D_CA, 303_D_C, 304_K_N)
     [  13.40836, 122.67608,  14.72403, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 303, "CNCA" ], // (303_D_C, 304_K_N, 304_K_CA)
     [  14.72403, 112.55464,  15.31465, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 304, "NCAC" ], // (304_K_N, 304_K_CA, 304_K_C)
     [  15.31465, 119.97917,  12.31634, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 304, "CACO" ], // (304_K_CA, 304_K_C, 304_K_O)
     [  15.31465, 118.29774,  13.34255, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 304, "CACN" ], // (304_K_CA, 304_K_C, 305_E_N)
     [  13.34255, 122.26469,  14.65759, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 304, "CNCA" ], // (304_K_C, 305_E_N, 305_E_CA)
     [  14.65759, 111.87771,  15.36852, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 305, "NCAC" ], // (305_E_N, 305_E_CA, 305_E_C)
     [  15.36852, 118.61158,  12.29059, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 305, "CACO" ], // (305_E_CA, 305_E_C, 305_E_O)
     [  15.36852, 119.68433,  13.39321, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 305, "CACN" ], // (305_E_CA, 305_E_C, 306_L_N)
     [  13.39321, 124.44527,  14.67871, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 305, "CNCA" ], // (305_E_C, 306_L_N, 306_L_CA)
     [  14.67871, 112.12192,  15.35606, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 306, "NCAC" ], // (306_L_N, 306_L_CA, 306_L_C)
     [  15.35606, 119.69695,  12.31834, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 306, "CACO" ], // (306_L_CA, 306_L_C, 306_L_O)
     [  15.35606, 118.80787,  13.37939, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 306, "CACN" ], // (306_L_CA, 306_L_C, 307_V_N)
     [  13.37939, 122.83788,  14.79602, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 306, "CNCA" ], // (306_L_C, 307_V_N, 307_V_CA)
     [  14.79602, 111.10054,  15.36130, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 307, "NCAC" ], // (307_V_N, 307_V_CA, 307_V_C)
     [  15.36130, 121.47779,  12.39425, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 307, "CACO" ], // (307_V_CA, 307_V_C, 307_V_O)
     [  15.36130, 116.94963,  13.33922, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 307, "CACN" ], // (307_V_CA, 307_V_C, 308_H_N)
     [  13.33922, 123.13509,  14.60218, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 307, "CNCA" ], // (307_V_C, 308_H_N, 308_H_CA)
     [  14.60218, 113.16380,  15.39305, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 308, "NCAC" ], // (308_H_N, 308_H_CA, 308_H_C)
     [  15.39305, 118.77348,  12.29774, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 308, "CACO" ], // (308_H_CA, 308_H_C, 308_H_O)
     [  15.39305, 119.53048,  13.36435, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 308, "CACN" ], // (308_H_CA, 308_H_C, 309_M_N)
     [  13.36435, 125.04492,  14.67774, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 308, "CNCA" ], // (308_H_C, 309_M_N, 309_M_CA)
     [  14.67774, 111.92334,  15.32874, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 309, "NCAC" ], // (309_M_N, 309_M_CA, 309_M_C)
     [  15.32874, 119.56633,  12.29080, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 309, "CACO" ], // (309_M_CA, 309_M_C, 309_M_O)
     [  15.32874, 118.66296,  13.36175, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 309, "CACN" ], // (309_M_CA, 309_M_C, 310_I_N)
     [  13.36175, 123.23849,  14.73835, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 309, "CNCA" ], // (309_M_C, 310_I_N, 310_I_CA)
     [  14.73835, 110.08208,  15.36316, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 310, "NCAC" ], // (310_I_N, 310_I_CA, 310_I_C)
     [  15.36316, 121.13062,  12.33091, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 310, "CACO" ], // (310_I_CA, 310_I_C, 310_I_O)
     [  15.36316, 116.79153,  13.35056, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 310, "CACN" ], // (310_I_CA, 310_I_C, 311_S_N)
     [  13.35056, 123.66945,  14.65619, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 310, "CNCA" ], // (310_I_C, 311_S_N, 311_S_CA)
     [  14.65619, 112.26897,  15.35607, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 311, "NCAC" ], // (311_S_N, 311_S_CA, 311_S_C)
     [  15.35607, 119.44996,  12.26607, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 311, "CACO" ], // (311_S_CA, 311_S_C, 311_S_O)
     [  15.35607, 118.31235,  13.36076, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 311, "CACN" ], // (311_S_CA, 311_S_C, 312_W_N)
     [  13.36076, 125.83173,  14.69039, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 311, "CNCA" ], // (311_S_C, 312_W_N, 312_W_CA)
     [  14.69039, 111.98373,  15.32523, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "W", 312, "NCAC" ], // (312_W_N, 312_W_CA, 312_W_C)
     [  15.32523, 119.35087,  12.27701, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "W", 312, "CACO" ], // (312_W_CA, 312_W_C, 312_W_O)
     [  15.32523, 118.44486,  13.36918, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "W", 312, "CACN" ], // (312_W_CA, 312_W_C, 313_A_N)
     [  13.36918, 123.73137,  14.62909, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "W", 312, "CNCA" ], // (312_W_C, 313_A_N, 313_A_CA)
     [  14.62909, 112.06094,  15.40552, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 313, "NCAC" ], // (313_A_N, 313_A_CA, 313_A_C)
     [  15.40552, 120.32436,  12.26733, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 313, "CACO" ], // (313_A_CA, 313_A_C, 313_A_O)
     [  15.40552, 117.85484,  13.41728, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 313, "CACN" ], // (313_A_CA, 313_A_C, 314_K_N)
     [  13.41728, 124.07450,  14.75580, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 313, "CNCA" ], // (313_A_C, 314_K_N, 314_K_CA)
     [  14.75580, 114.75732,  15.36012, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 314, "NCAC" ], // (314_K_N, 314_K_CA, 314_K_C)
     [  15.36012, 118.39405,  12.27463, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 314, "CACO" ], // (314_K_CA, 314_K_C, 314_K_O)
     [  15.36012, 119.25039,  13.36961, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 314, "CACN" ], // (314_K_CA, 314_K_C, 315_K_N)
     [  13.36961, 123.85505,  14.73618, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 314, "CNCA" ], // (314_K_C, 315_K_N, 315_K_CA)
     [  14.73618, 113.00706,  15.33629, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 315, "NCAC" ], // (315_K_N, 315_K_CA, 315_K_C)
     [  15.33629, 118.76686,  12.27623, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 315, "CACO" ], // (315_K_CA, 315_K_C, 315_K_O)
     [  15.33629, 119.14390,  13.36695, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 315, "CACN" ], // (315_K_CA, 315_K_C, 316_I_N)
     [  13.36695, 122.15318,  14.66323, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 315, "CNCA" ], // (315_K_C, 316_I_N, 316_I_CA)
     [  14.66323, 108.39265,  15.37247, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 316, "NCAC" ], // (316_I_N, 316_I_CA, 316_I_C)
     [  15.37247, 118.96853,  12.35659, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 316, "CACO" ], // (316_I_CA, 316_I_C, 316_I_O)
     [  15.37247, 119.15542,  13.49650, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 316, "CACN" ], // (316_I_CA, 316_I_C, 317_P_N)
     [  13.49650, 122.43155,  14.58203, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 316, "CNCA" ], // (316_I_C, 317_P_N, 317_P_CA)
     [  14.58203, 112.06881,  15.31188, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 317, "NCAC" ], // (317_P_N, 317_P_CA, 317_P_C)
     [  15.31188, 120.72542,  12.31183, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 317, "CACO" ], // (317_P_CA, 317_P_C, 317_P_O)
     [  15.31188, 117.09370,  13.33993, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 317, "CACN" ], // (317_P_CA, 317_P_C, 318_G_N)
     [  13.33993, 123.37016,  14.62105, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 317, "CNCA" ], // (317_P_C, 318_G_N, 318_G_CA)
     [  14.62105, 115.59754,  15.31508, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 318, "NCAC" ], // (318_G_N, 318_G_CA, 318_G_C)
     [  15.31508, 118.78973,  12.30409, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 318, "CACO" ], // (318_G_CA, 318_G_C, 318_G_O)
     [  15.31508, 119.62743,  13.39425, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 318, "CACN" ], // (318_G_CA, 318_G_C, 319_F_N)
     [  13.39425, 122.26027,  14.67797, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 318, "CNCA" ], // (318_G_C, 319_F_N, 319_F_CA)
     [  14.67797, 112.02544,  15.35944, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 319, "NCAC" ], // (319_F_N, 319_F_CA, 319_F_C)
     [  15.35944, 120.29595,  12.30436, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 319, "CACO" ], // (319_F_CA, 319_F_C, 319_F_O)
     [  15.35944, 118.03125,  13.38769, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 319, "CACN" ], // (319_F_CA, 319_F_C, 320_V_N)
     [  13.38769, 123.17769,  14.76395, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 319, "CNCA" ], // (319_F_C, 320_V_N, 320_V_CA)
     [  14.76395, 113.16019,  15.34779, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 320, "NCAC" ], // (320_V_N, 320_V_CA, 320_V_C)
     [  15.34779, 120.15097,  12.30498, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 320, "CACO" ], // (320_V_CA, 320_V_C, 320_V_O)
     [  15.34779, 117.42606,  13.34050, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 320, "CACN" ], // (320_V_CA, 320_V_C, 321_E_N)
     [  13.34050, 126.91611,  14.64428, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 320, "CNCA" ], // (320_V_C, 321_E_N, 321_E_CA)
     [  14.64428, 113.60935,  15.43554, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 321, "NCAC" ], // (321_E_N, 321_E_CA, 321_E_C)
     [  15.43554, 119.61894,  12.37450, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 321, "CACO" ], // (321_E_CA, 321_E_C, 321_E_O)
     [  15.43554, 118.39348,  13.39296, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 321, "CACN" ], // (321_E_CA, 321_E_C, 322_L_N)
     [  13.39296, 123.82972,  14.67576, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 321, "CNCA" ], // (321_E_C, 322_L_N, 322_L_CA)
     [  14.67576, 110.46553,  15.37457, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 322, "NCAC" ], // (322_L_N, 322_L_CA, 322_L_C)
     [  15.37457, 120.77831,  12.33608, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 322, "CACO" ], // (322_L_CA, 322_L_C, 322_L_O)
     [  15.37457, 118.36316,  13.37188, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 322, "CACN" ], // (322_L_CA, 322_L_C, 323_S_N)
     [  13.37188, 120.33968,  14.68799, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 322, "CNCA" ], // (322_L_C, 323_S_N, 323_S_CA)
     [  14.68799, 110.86289,  15.35775, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 323, "NCAC" ], // (323_S_N, 323_S_CA, 323_S_C)
     [  15.35775, 120.37724,  12.30077, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 323, "CACO" ], // (323_S_CA, 323_S_C, 323_S_O)
     [  15.35775, 118.33581,  13.39067, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 323, "CACN" ], // (323_S_CA, 323_S_C, 324_L_N)
     [  13.39067, 121.81330,  14.73405, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 323, "CNCA" ], // (323_S_C, 324_L_N, 324_L_CA)
     [  14.73405, 112.09374,  15.37928, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 324, "NCAC" ], // (324_L_N, 324_L_CA, 324_L_C)
     [  15.37928, 119.88774,  12.27217, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 324, "CACO" ], // (324_L_CA, 324_L_C, 324_L_O)
     [  15.37928, 117.29342,  13.36738, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 324, "CACN" ], // (324_L_CA, 324_L_C, 325_F_N)
     [  13.36738, 123.58703,  14.61425, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 324, "CNCA" ], // (324_L_C, 325_F_N, 325_F_CA)
     [  14.61425, 111.61954,  15.32585, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 325, "NCAC" ], // (325_F_N, 325_F_CA, 325_F_C)
     [  15.32585, 119.82030,  12.27880, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 325, "CACO" ], // (325_F_CA, 325_F_C, 325_F_O)
     [  15.32585, 117.86184,  13.35200, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 325, "CACN" ], // (325_F_CA, 325_F_C, 326_D_N)
     [  13.35200, 122.66104,  14.64767, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 325, "CNCA" ], // (325_F_C, 326_D_N, 326_D_CA)
     [  14.64767, 111.17518,  15.32872, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 326, "NCAC" ], // (326_D_N, 326_D_CA, 326_D_C)
     [  15.32872, 119.26595,  12.26214, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 326, "CACO" ], // (326_D_CA, 326_D_C, 326_D_O)
     [  15.32872, 118.61098,  13.38569, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 326, "CACN" ], // (326_D_CA, 326_D_C, 327_Q_N)
     [  13.38569, 123.62224,  14.70222, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 326, "CNCA" ], // (326_D_C, 327_Q_N, 327_Q_CA)
     [  14.70222, 111.71212,  15.36594, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 327, "NCAC" ], // (327_Q_N, 327_Q_CA, 327_Q_C)
     [  15.36594, 119.84744,  12.32257, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 327, "CACO" ], // (327_Q_CA, 327_Q_C, 327_Q_O)
     [  15.36594, 117.80433,  13.39919, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 327, "CACN" ], // (327_Q_CA, 327_Q_C, 328_V_N)
     [  13.39919, 123.36595,  14.73356, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 327, "CNCA" ], // (327_Q_C, 328_V_N, 328_V_CA)
     [  14.73356, 110.46398,  15.27957, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 328, "NCAC" ], // (328_V_N, 328_V_CA, 328_V_C)
     [  15.27957, 119.81464,  12.28092, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 328, "CACO" ], // (328_V_CA, 328_V_C, 328_V_O)
     [  15.27957, 118.41308,  13.38546, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 328, "CACN" ], // (328_V_CA, 328_V_C, 329_R_N)
     [  13.38546, 123.61167,  14.68264, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 328, "CNCA" ], // (328_V_C, 329_R_N, 329_R_CA)
     [  14.68264, 113.10293,  15.33183, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 329, "NCAC" ], // (329_R_N, 329_R_CA, 329_R_C)
     [  15.33183, 117.98995,  12.36930, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 329, "CACO" ], // (329_R_CA, 329_R_C, 329_R_O)
     [  15.33183, 117.93914,  13.38074, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 329, "CACN" ], // (329_R_CA, 329_R_C, 330_L_N)
     [  13.38074, 125.37429,  14.66545, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 329, "CNCA" ], // (329_R_C, 330_L_N, 330_L_CA)
     [  14.66545, 110.51386,  15.30599, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 330, "NCAC" ], // (330_L_N, 330_L_CA, 330_L_C)
     [  15.30599, 119.05784,  12.30155, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 330, "CACO" ], // (330_L_CA, 330_L_C, 330_L_O)
     [  15.30599, 118.06102,  13.36881, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 330, "CACN" ], // (330_L_CA, 330_L_C, 331_L_N)
     [  13.36881, 124.09278,  14.69123, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 330, "CNCA" ], // (330_L_C, 331_L_N, 331_L_CA)
     [  14.69123, 111.23439,  15.38985, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 331, "NCAC" ], // (331_L_N, 331_L_CA, 331_L_C)
     [  15.38985, 120.00236,  12.29706, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 331, "CACO" ], // (331_L_CA, 331_L_C, 331_L_O)
     [  15.38985, 117.42575,  13.36005, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 331, "CACN" ], // (331_L_CA, 331_L_C, 332_E_N)
     [  13.36005, 125.30561,  14.70560, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 331, "CNCA" ], // (331_L_C, 332_E_N, 332_E_CA)
     [  14.70560, 113.63587,  15.34875, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 332, "NCAC" ], // (332_E_N, 332_E_CA, 332_E_C)
     [  15.34875, 119.17862,  12.21857, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 332, "CACO" ], // (332_E_CA, 332_E_C, 332_E_O)
     [  15.34875, 119.73624,  13.38765, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 332, "CACN" ], // (332_E_CA, 332_E_C, 333_S_N)
     [  13.38765, 122.51900,  14.70412, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 332, "CNCA" ], // (332_E_C, 333_S_N, 333_S_CA)
     [  14.70412, 114.59463,  15.26149, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 333, "NCAC" ], // (333_S_N, 333_S_CA, 333_S_C)
     [  15.26149, 118.82258,  12.28945, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 333, "CACO" ], // (333_S_CA, 333_S_C, 333_S_O)
     [  15.26149, 117.19472,  13.35277, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 333, "CACN" ], // (333_S_CA, 333_S_C, 334_C_N)
     [  13.35277, 123.65402,  14.63343, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 333, "CNCA" ], // (333_S_C, 334_C_N, 334_C_CA)
     [  14.63343, 113.77886,  15.33389, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "C", 334, "NCAC" ], // (334_C_N, 334_C_CA, 334_C_C)
     [  15.33389, 119.66318,  12.28658, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "C", 334, "CACO" ], // (334_C_CA, 334_C_C, 334_C_O)
     [  15.33389, 119.29566,  13.35565, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "C", 334, "CACN" ], // (334_C_CA, 334_C_C, 335_W_N)
     [  13.35565, 121.79934,  14.63484, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "C", 334, "CNCA" ], // (334_C_C, 335_W_N, 335_W_CA)
     [  14.63484, 113.46348,  15.28933, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "W", 335, "NCAC" ], // (335_W_N, 335_W_CA, 335_W_C)
     [  15.28933, 119.76282,  12.27524, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "W", 335, "CACO" ], // (335_W_CA, 335_W_C, 335_W_O)
     [  15.28933, 117.62532,  13.32093, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "W", 335, "CACN" ], // (335_W_CA, 335_W_C, 336_M_N)
     [  13.32093, 122.96038,  14.62254, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "W", 335, "CNCA" ], // (335_W_C, 336_M_N, 336_M_CA)
     [  14.62254, 111.83948,  15.39837, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 336, "NCAC" ], // (336_M_N, 336_M_CA, 336_M_C)
     [  15.39837, 119.69730,  12.28114, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 336, "CACO" ], // (336_M_CA, 336_M_C, 336_M_O)
     [  15.39837, 118.75868,  13.37393, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 336, "CACN" ], // (336_M_CA, 336_M_C, 337_E_N)
     [  13.37393, 123.24529,  14.68432, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 336, "CNCA" ], // (336_M_C, 337_E_N, 337_E_CA)
     [  14.68432, 111.84955,  15.31968, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 337, "NCAC" ], // (337_E_N, 337_E_CA, 337_E_C)
     [  15.31968, 119.67809,  12.28719, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 337, "CACO" ], // (337_E_CA, 337_E_C, 337_E_O)
     [  15.31968, 118.39319,  13.37402, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 337, "CACN" ], // (337_E_CA, 337_E_C, 338_V_N)
     [  13.37402, 123.55297,  14.71936, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 337, "CNCA" ], // (337_E_C, 338_V_N, 338_V_CA)
     [  14.71936, 111.27331,  15.31859, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 338, "NCAC" ], // (338_V_N, 338_V_CA, 338_V_C)
     [  15.31859, 118.93573,  12.34866, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 338, "CACO" ], // (338_V_CA, 338_V_C, 338_V_O)
     [  15.31859, 118.65807,  13.39505, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 338, "CACN" ], // (338_V_CA, 338_V_C, 339_L_N)
     [  13.39505, 123.05415,  14.67968, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 338, "CNCA" ], // (338_V_C, 339_L_N, 339_L_CA)
     [  14.67968, 111.52757,  15.35234, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 339, "NCAC" ], // (339_L_N, 339_L_CA, 339_L_C)
     [  15.35234, 119.99066,  12.30153, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 339, "CACO" ], // (339_L_CA, 339_L_C, 339_L_O)
     [  15.35234, 118.39513,  13.38655, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 339, "CACN" ], // (339_L_CA, 339_L_C, 340_M_N)
     [  13.38655, 122.77739,  14.67987, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 339, "CNCA" ], // (339_L_C, 340_M_N, 340_M_CA)
     [  14.67987, 112.50208,  15.33382, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 340, "NCAC" ], // (340_M_N, 340_M_CA, 340_M_C)
     [  15.33382, 120.01572,  12.28849, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 340, "CACO" ], // (340_M_CA, 340_M_C, 340_M_O)
     [  15.33382, 118.43069,  13.35484, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 340, "CACN" ], // (340_M_CA, 340_M_C, 341_M_N)
     [  13.35484, 123.57963,  14.64637, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 340, "CNCA" ], // (340_M_C, 341_M_N, 341_M_CA)
     [  14.64637, 111.04658,  15.38467, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 341, "NCAC" ], // (341_M_N, 341_M_CA, 341_M_C)
     [  15.38467, 119.90500,  12.32656, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 341, "CACO" ], // (341_M_CA, 341_M_C, 341_M_O)
     [  15.38467, 117.77789,  13.39334, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 341, "CACN" ], // (341_M_CA, 341_M_C, 342_G_N)
     [  13.39334, 123.93171,  14.65216, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 341, "CNCA" ], // (341_M_C, 342_G_N, 342_G_CA)
     [  14.65216, 113.34560,  15.32364, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 342, "NCAC" ], // (342_G_N, 342_G_CA, 342_G_C)
     [  15.32364, 119.81685,  12.30941, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 342, "CACO" ], // (342_G_CA, 342_G_C, 342_G_O)
     [  15.32364, 117.95764,  13.34766, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 342, "CACN" ], // (342_G_CA, 342_G_C, 343_L_N)
     [  13.34766, 124.61408,  14.65115, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 342, "CNCA" ], // (342_G_C, 343_L_N, 343_L_CA)
     [  14.65115, 112.27580,  15.34386, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 343, "NCAC" ], // (343_L_N, 343_L_CA, 343_L_C)
     [  15.34386, 120.33514,  12.31431, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 343, "CACO" ], // (343_L_CA, 343_L_C, 343_L_O)
     [  15.34386, 117.74718,  13.34326, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 343, "CACN" ], // (343_L_CA, 343_L_C, 344_M_N)
     [  13.34326, 122.79904,  14.67116, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 343, "CNCA" ], // (343_L_C, 344_M_N, 344_M_CA)
     [  14.67116, 111.98635,  15.32069, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 344, "NCAC" ], // (344_M_N, 344_M_CA, 344_M_C)
     [  15.32069, 119.82949,  12.27895, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 344, "CACO" ], // (344_M_CA, 344_M_C, 344_M_O)
     [  15.32069, 118.12854,  13.35426, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 344, "CACN" ], // (344_M_CA, 344_M_C, 345_W_N)
     [  13.35426, 122.81382,  14.64282, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 344, "CNCA" ], // (344_M_C, 345_W_N, 345_W_CA)
     [  14.64282, 112.47997,  15.32030, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "W", 345, "NCAC" ], // (345_W_N, 345_W_CA, 345_W_C)
     [  15.32030, 119.50312,  12.27350, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "W", 345, "CACO" ], // (345_W_CA, 345_W_C, 345_W_O)
     [  15.32030, 118.39977,  13.42316, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "W", 345, "CACN" ], // (345_W_CA, 345_W_C, 346_R_N)
     [  13.42316, 123.22641,  14.79993, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "W", 345, "CNCA" ], // (345_W_C, 346_R_N, 346_R_CA)
     [  14.79993, 113.35715,  15.32585, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 346, "NCAC" ], // (346_R_N, 346_R_CA, 346_R_C)
     [  15.32585, 119.90247,  12.25172, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 346, "CACO" ], // (346_R_CA, 346_R_C, 346_R_O)
     [  15.32585, 117.52778,  13.32667, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 346, "CACN" ], // (346_R_CA, 346_R_C, 347_S_N)
     [  13.32667, 124.64248,  14.69734, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 346, "CNCA" ], // (346_R_C, 347_S_N, 347_S_CA)
     [  14.69734, 113.98262,  15.32045, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 347, "NCAC" ], // (347_S_N, 347_S_CA, 347_S_C)
     [  15.32045, 119.44595,  12.33173, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 347, "CACO" ], // (347_S_CA, 347_S_C, 347_S_O)
     [  15.32045, 118.94203,  13.40397, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 347, "CACN" ], // (347_S_CA, 347_S_C, 348_I_N)
     [  13.40397, 122.76764,  14.68861, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 347, "CNCA" ], // (347_S_C, 348_I_N, 348_I_CA)
     [  14.68861, 112.09981,  15.36743, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 348, "NCAC" ], // (348_I_N, 348_I_CA, 348_I_C)
     [  15.36743, 119.21049,  12.34222, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 348, "CACO" ], // (348_I_CA, 348_I_C, 348_I_O)
     [  15.36743, 117.67806,  13.34543, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 348, "CACN" ], // (348_I_CA, 348_I_C, 349_D_N)
     [  13.34543, 125.54364,  14.71489, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 348, "CNCA" ], // (348_I_C, 349_D_N, 349_D_CA)
     [  14.71489, 114.40361,  15.43895, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 349, "NCAC" ], // (349_D_N, 349_D_CA, 349_D_C)
     [  15.43895, 120.71990,  12.31456, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 349, "CACO" ], // (349_D_CA, 349_D_C, 349_D_O)
     [  15.43895, 118.32727,  13.39080, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 349, "CACN" ], // (349_D_CA, 349_D_C, 350_H_N)
     [  13.39080, 124.27209,  14.67602, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 349, "CNCA" ], // (349_D_C, 350_H_N, 350_H_CA)
     [  14.67602, 112.52133,  15.45065, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 350, "NCAC" ], // (350_H_N, 350_H_CA, 350_H_C)
     [  15.45065, 120.69859,  12.32014, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 350, "CACO" ], // (350_H_CA, 350_H_C, 350_H_O)
     [  15.45065, 119.38525,  13.50273, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 350, "CACN" ], // (350_H_CA, 350_H_C, 351_P_N)
     [  13.50273, 122.76916,  14.64269, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 350, "CNCA" ], // (350_H_C, 351_P_N, 351_P_CA)
     [  14.64269, 113.08334,  15.37399, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 351, "NCAC" ], // (351_P_N, 351_P_CA, 351_P_C)
     [  15.37399, 120.97751,  12.28838, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 351, "CACO" ], // (351_P_CA, 351_P_C, 351_P_O)
     [  15.37399, 116.59314,  13.36773, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 351, "CACN" ], // (351_P_CA, 351_P_C, 352_G_N)
     [  13.36773, 124.15239,  14.68422, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 351, "CNCA" ], // (351_P_C, 352_G_N, 352_G_CA)
     [  14.68422, 115.81111,  15.36022, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 352, "NCAC" ], // (352_G_N, 352_G_CA, 352_G_C)
     [  15.36022, 118.51663,  12.28119, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 352, "CACO" ], // (352_G_CA, 352_G_C, 352_G_O)
     [  15.36022, 119.16817,  13.40766, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 352, "CACN" ], // (352_G_CA, 352_G_C, 353_K_N)
     [  13.40766, 124.41328,  14.72473, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 352, "CNCA" ], // (352_G_C, 353_K_N, 353_K_CA)
     [  14.72473, 111.53069,  15.26634, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 353, "NCAC" ], // (353_K_N, 353_K_CA, 353_K_C)
     [  15.26634, 120.67199,  12.28627, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 353, "CACO" ], // (353_K_CA, 353_K_C, 353_K_O)
     [  15.26634, 116.28354,  13.29304, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 353, "CACN" ], // (353_K_CA, 353_K_C, 354_L_N)
     [  13.29304, 124.85238,  14.58834, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 353, "CNCA" ], // (353_K_C, 354_L_N, 354_L_CA)
     [  14.58834, 109.94516,  15.35630, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 354, "NCAC" ], // (354_L_N, 354_L_CA, 354_L_C)
     [  15.35630, 121.58552,  12.31981, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 354, "CACO" ], // (354_L_CA, 354_L_C, 354_L_O)
     [  15.35630, 115.92807,  13.32876, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 354, "CACN" ], // (354_L_CA, 354_L_C, 355_I_N)
     [  13.32876, 125.79491,  14.66903, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 354, "CNCA" ], // (354_L_C, 355_I_N, 355_I_CA)
     [  14.66903, 108.66421,  15.33326, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 355, "NCAC" ], // (355_I_N, 355_I_CA, 355_I_C)
     [  15.33326, 122.06837,  12.37077, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 355, "CACO" ], // (355_I_CA, 355_I_C, 355_I_O)
     [  15.33326, 115.67634,  13.33818, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 355, "CACN" ], // (355_I_CA, 355_I_C, 356_F_N)
     [  13.33818, 123.06733,  14.63892, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 355, "CNCA" ], // (355_I_C, 356_F_N, 356_F_CA)
     [  14.63892, 110.84273,  15.27574, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 356, "NCAC" ], // (356_F_N, 356_F_CA, 356_F_C)
     [  15.27574, 119.51086,  12.27900, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 356, "CACO" ], // (356_F_CA, 356_F_C, 356_F_O)
     [  15.27574, 117.48107,  13.27774, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 356, "CACN" ], // (356_F_CA, 356_F_C, 357_A_N)
     [  13.27774, 124.76411,  14.59171, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 356, "CNCA" ], // (356_F_C, 357_A_N, 357_A_CA)
     [  14.59171, 109.61677,  15.42569, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 357, "NCAC" ], // (357_A_N, 357_A_CA, 357_A_C)
     [  15.42569, 118.72385,  12.32012, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 357, "CACO" ], // (357_A_CA, 357_A_C, 357_A_O)
     [  15.42569, 120.33739,  13.48211, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 357, "CACN" ], // (357_A_CA, 357_A_C, 358_P_N)
     [  13.48211, 122.44405,  14.67736, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 357, "CNCA" ], // (357_A_C, 358_P_N, 358_P_CA)
     [  14.67736, 113.65855,  15.38537, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 358, "NCAC" ], // (358_P_N, 358_P_CA, 358_P_C)
     [  15.38537, 120.12038,  12.27314, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 358, "CACO" ], // (358_P_CA, 358_P_C, 358_P_O)
     [  15.38537, 117.57631,  13.37589, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 358, "CACN" ], // (358_P_CA, 358_P_C, 359_D_N)
     [  13.37589, 125.99466,  14.64141, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 358, "CNCA" ], // (358_P_C, 359_D_N, 359_D_CA)
     [  14.64141, 114.10915,  15.43283, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 359, "NCAC" ], // (359_D_N, 359_D_CA, 359_D_C)
     [  15.43283, 119.27750,  12.31162, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 359, "CACO" ], // (359_D_CA, 359_D_C, 359_D_O)
     [  15.43283, 119.36390,  13.39041, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 359, "CACN" ], // (359_D_CA, 359_D_C, 360_L_N)
     [  13.39041, 125.01652,  14.70050, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 359, "CNCA" ], // (359_D_C, 360_L_N, 360_L_CA)
     [  14.70050, 110.76477,  15.34615, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 360, "NCAC" ], // (360_L_N, 360_L_CA, 360_L_C)
     [  15.34615, 121.93347,  12.32034, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 360, "CACO" ], // (360_L_CA, 360_L_C, 360_L_O)
     [  15.34615, 114.89256,  13.30564, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 360, "CACN" ], // (360_L_CA, 360_L_C, 361_V_N)
     [  13.30564, 124.57358,  14.65951, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 360, "CNCA" ], // (360_L_C, 361_V_N, 361_V_CA)
     [  14.65951, 108.26024,  15.17881, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 361, "NCAC" ], // (361_V_N, 361_V_CA, 361_V_C)
     [  15.17881, 118.51111,  12.29922, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 361, "CACO" ], // (361_V_CA, 361_V_C, 361_V_O)
     [  15.17881, 116.73172,  13.32381, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 361, "CACN" ], // (361_V_CA, 361_V_C, 362_L_N)
     [  13.32381, 124.36907,  14.59498, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 361, "CNCA" ], // (361_V_C, 362_L_N, 362_L_CA)
     [  14.59498, 110.45593,  15.27202, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 362, "NCAC" ], // (362_L_N, 362_L_CA, 362_L_C)
     [  15.27202, 120.40348,  12.34769, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 362, "CACO" ], // (362_L_CA, 362_L_C, 362_L_O)
     [  15.27202, 117.16199,  13.26830, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 362, "CACN" ], // (362_L_CA, 362_L_C, 363_D_N)
     [  13.26830, 123.35366,  14.62317, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 362, "CNCA" ], // (362_L_C, 363_D_N, 363_D_CA)
     [  14.62317, 112.61388,  15.41724, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 363, "NCAC" ], // (363_D_N, 363_D_CA, 363_D_C)
     [  15.41724, 122.42955,  12.35810, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 363, "CACO" ], // (363_D_CA, 363_D_C, 363_D_O)
     [  15.41724, 115.80968,  13.37916, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 363, "CACN" ], // (363_D_CA, 363_D_C, 364_R_N)
     [  13.37916, 123.37229,  14.67991, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 363, "CNCA" ], // (363_D_C, 364_R_N, 364_R_CA)
     [  14.67991, 112.94104,  15.31716, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 364, "NCAC" ], // (364_R_N, 364_R_CA, 364_R_C)
     [  15.31716, 118.78668,  12.28972, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 364, "CACO" ], // (364_R_CA, 364_R_C, 364_R_O)
     [  15.31716, 118.86140,  13.32270, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 364, "CACN" ], // (364_R_CA, 364_R_C, 365_D_N)
     [  13.32270, 123.19030,  14.62045, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 364, "CNCA" ], // (364_R_C, 365_D_N, 365_D_CA)
     [  14.62045, 110.82127,  15.31320, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 365, "NCAC" ], // (365_D_N, 365_D_CA, 365_D_C)
     [  15.31320, 119.45428,  12.30087, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 365, "CACO" ], // (365_D_CA, 365_D_C, 365_D_O)
     [  15.31320, 118.12096,  13.35372, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 365, "CACN" ], // (365_D_CA, 365_D_C, 366_E_N)
     [  13.35372, 123.58560,  14.71246, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 365, "CNCA" ], // (365_D_C, 366_E_N, 366_E_CA)
     [  14.71246, 112.14953,  15.40029, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 366, "NCAC" ], // (366_E_N, 366_E_CA, 366_E_C)
     [  15.40029, 120.51247,  12.34166, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 366, "CACO" ], // (366_E_CA, 366_E_C, 366_E_O)
     [  15.40029, 117.72778,  13.37975, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 366, "CACN" ], // (366_E_CA, 366_E_C, 367_G_N)
     [  13.37975, 123.89394,  14.66144, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 366, "CNCA" ], // (366_E_C, 367_G_N, 367_G_CA)
     [  14.66144, 115.07258,  15.37665, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 367, "NCAC" ], // (367_G_N, 367_G_CA, 367_G_C)
     [  15.37665, 119.23461,  12.28472, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 367, "CACO" ], // (367_G_CA, 367_G_C, 367_G_O)
     [  15.37665, 119.96574,  13.42689, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 367, "CACN" ], // (367_G_CA, 367_G_C, 368_K_N)
     [  13.42689, 123.38066,  14.81774, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 367, "CNCA" ], // (367_G_C, 368_K_N, 368_K_CA)
     [  14.81774, 114.35326,  15.30629, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 368, "NCAC" ], // (368_K_N, 368_K_CA, 368_K_C)
     [  15.30629, 118.94767,  12.26345, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 368, "CACO" ], // (368_K_CA, 368_K_C, 368_K_O)
     [  15.30629, 118.57577,  13.33286, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 368, "CACN" ], // (368_K_CA, 368_K_C, 369_C_N)
     [  13.33286, 123.07863,  14.61913, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 368, "CNCA" ], // (368_K_C, 369_C_N, 369_C_CA)
     [  14.61913, 114.53743,  15.41019, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "C", 369, "NCAC" ], // (369_C_N, 369_C_CA, 369_C_C)
     [  15.41019, 117.68008,  12.26209, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "C", 369, "CACO" ], // (369_C_CA, 369_C_C, 369_C_O)
     [  15.41019, 119.79924,  13.40646, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "C", 369, "CACN" ], // (369_C_CA, 369_C_C, 370_V_N)
     [  13.40646, 124.58185,  14.75682, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "C", 369, "CNCA" ], // (369_C_C, 370_V_N, 370_V_CA)
     [  14.75682, 110.39909,  15.30324, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 370, "NCAC" ], // (370_V_N, 370_V_CA, 370_V_C)
     [  15.30324, 120.68081,  12.35340, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 370, "CACO" ], // (370_V_CA, 370_V_C, 370_V_O)
     [  15.30324, 117.66633,  13.34005, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 370, "CACN" ], // (370_V_CA, 370_V_C, 371_E_N)
     [  13.34005, 122.33309,  14.66108, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 370, "CNCA" ], // (370_V_C, 371_E_N, 371_E_CA)
     [  14.66108, 111.79628,  15.42238, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 371, "NCAC" ], // (371_E_N, 371_E_CA, 371_E_C)
     [  15.42238, 122.16282,  12.30838, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 371, "CACO" ], // (371_E_CA, 371_E_C, 371_E_O)
     [  15.42238, 116.94610,  13.36565, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 371, "CACN" ], // (371_E_CA, 371_E_C, 372_G_N)
     [  13.36565, 123.81773,  14.68088, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 371, "CNCA" ], // (371_E_C, 372_G_N, 372_G_CA)
     [  14.68088, 116.04703,  15.35536, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 372, "NCAC" ], // (372_G_N, 372_G_CA, 372_G_C)
     [  15.35536, 118.66773,  12.28877, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 372, "CACO" ], // (372_G_CA, 372_G_C, 372_G_O)
     [  15.35536, 120.20673,  13.41404, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 372, "CACN" ], // (372_G_CA, 372_G_C, 373_I_N)
     [  13.41404, 121.77751,  14.77761, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 372, "CNCA" ], // (372_G_C, 373_I_N, 373_I_CA)
     [  14.77761, 113.02267,  15.39794, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 373, "NCAC" ], // (373_I_N, 373_I_CA, 373_I_C)
     [  15.39794, 120.32041,  12.30122, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 373, "CACO" ], // (373_I_CA, 373_I_C, 373_I_O)
     [  15.39794, 118.00931,  13.39340, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 373, "CACN" ], // (373_I_CA, 373_I_C, 374_L_N)
     [  13.39340, 123.85510,  14.68364, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 373, "CNCA" ], // (373_I_C, 374_L_N, 374_L_CA)
     [  14.68364, 113.33955,  15.35307, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 374, "NCAC" ], // (374_L_N, 374_L_CA, 374_L_C)
     [  15.35307, 118.48865,  12.26713, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 374, "CACO" ], // (374_L_CA, 374_L_C, 374_L_O)
     [  15.35307, 119.16236,  13.39363, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 374, "CACN" ], // (374_L_CA, 374_L_C, 375_E_N)
     [  13.39363, 122.38111,  14.73186, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 374, "CNCA" ], // (374_L_C, 375_E_N, 375_E_CA)
     [  14.73186, 110.26742,  15.37640, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 375, "NCAC" ], // (375_E_N, 375_E_CA, 375_E_C)
     [  15.37640, 121.00151,  12.29709, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 375, "CACO" ], // (375_E_CA, 375_E_C, 375_E_O)
     [  15.37640, 117.06172,  13.37135, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 375, "CACN" ], // (375_E_CA, 375_E_C, 376_I_N)
     [  13.37135, 124.15421,  14.69640, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 375, "CNCA" ], // (375_E_C, 376_I_N, 376_I_CA)
     [  14.69640, 110.45145,  15.34579, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 376, "NCAC" ], // (376_I_N, 376_I_CA, 376_I_C)
     [  15.34579, 120.13917,  12.30077, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 376, "CACO" ], // (376_I_CA, 376_I_C, 376_I_O)
     [  15.34579, 117.78254,  13.35708, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 376, "CACN" ], // (376_I_CA, 376_I_C, 377_F_N)
     [  13.35708, 123.73254,  14.64051, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 376, "CNCA" ], // (376_I_C, 377_F_N, 377_F_CA)
     [  14.64051, 111.49004,  15.37247, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 377, "NCAC" ], // (377_F_N, 377_F_CA, 377_F_C)
     [  15.37247, 120.25345,  12.31906, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 377, "CACO" ], // (377_F_CA, 377_F_C, 377_F_O)
     [  15.37247, 117.57075,  13.35145, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 377, "CACN" ], // (377_F_CA, 377_F_C, 378_D_N)
     [  13.35145, 123.71059,  14.64521, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 377, "CNCA" ], // (377_F_C, 378_D_N, 378_D_CA)
     [  14.64521, 111.77425,  15.42520, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 378, "NCAC" ], // (378_D_N, 378_D_CA, 378_D_C)
     [  15.42520, 119.77225,  12.34099, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 378, "CACO" ], // (378_D_CA, 378_D_C, 378_D_O)
     [  15.42520, 118.92273,  13.40894, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 378, "CACN" ], // (378_D_CA, 378_D_C, 379_M_N)
     [  13.40894, 123.55902,  14.74350, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 378, "CNCA" ], // (378_D_C, 379_M_N, 379_M_CA)
     [  14.74350, 112.74793,  15.39646, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 379, "NCAC" ], // (379_M_N, 379_M_CA, 379_M_C)
     [  15.39646, 120.33957,  12.37874, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 379, "CACO" ], // (379_M_CA, 379_M_C, 379_M_O)
     [  15.39646, 117.64588,  13.37579, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 379, "CACN" ], // (379_M_CA, 379_M_C, 380_L_N)
     [  13.37579, 124.11278,  14.65677, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 379, "CNCA" ], // (379_M_C, 380_L_N, 380_L_CA)
     [  14.65677, 111.35903,  15.32945, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 380, "NCAC" ], // (380_L_N, 380_L_CA, 380_L_C)
     [  15.32945, 119.61728,  12.33343, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 380, "CACO" ], // (380_L_CA, 380_L_C, 380_L_O)
     [  15.32945, 118.70212,  13.38236, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 380, "CACN" ], // (380_L_CA, 380_L_C, 381_L_N)
     [  13.38236, 123.00325,  14.66909, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 380, "CNCA" ], // (380_L_C, 381_L_N, 381_L_CA)
     [  14.66909, 112.27779,  15.35405, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 381, "NCAC" ], // (381_L_N, 381_L_CA, 381_L_C)
     [  15.35405, 119.27090,  12.34052, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 381, "CACO" ], // (381_L_CA, 381_L_C, 381_L_O)
     [  15.35405, 118.39200,  13.36725, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 381, "CACN" ], // (381_L_CA, 381_L_C, 382_A_N)
     [  13.36725, 123.53512,  14.62238, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 381, "CNCA" ], // (381_L_C, 382_A_N, 382_A_CA)
     [  14.62238, 111.49627,  15.38844, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 382, "NCAC" ], // (382_A_N, 382_A_CA, 382_A_C)
     [  15.38844, 120.82512,  12.31879, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 382, "CACO" ], // (382_A_CA, 382_A_C, 382_A_O)
     [  15.38844, 117.58714,  13.35415, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 382, "CACN" ], // (382_A_CA, 382_A_C, 383_T_N)
     [  13.35415, 122.74426,  14.71643, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 382, "CNCA" ], // (382_A_C, 383_T_N, 383_T_CA)
     [  14.71643, 111.70764,  15.37099, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 383, "NCAC" ], // (383_T_N, 383_T_CA, 383_T_C)
     [  15.37099, 119.98329,  12.31618, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 383, "CACO" ], // (383_T_CA, 383_T_C, 383_T_O)
     [  15.37099, 118.57205,  13.38053, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 383, "CACN" ], // (383_T_CA, 383_T_C, 384_T_N)
     [  13.38053, 123.79205,  14.76700, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 383, "CNCA" ], // (383_T_C, 384_T_N, 384_T_CA)
     [  14.76700, 111.84165,  15.36228, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 384, "NCAC" ], // (384_T_N, 384_T_CA, 384_T_C)
     [  15.36228, 120.51232,  12.31842, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 384, "CACO" ], // (384_T_CA, 384_T_C, 384_T_O)
     [  15.36228, 118.80098,  13.40330, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 384, "CACN" ], // (384_T_CA, 384_T_C, 385_S_N)
     [  13.40330, 122.27089,  14.70431, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 384, "CNCA" ], // (384_T_C, 385_S_N, 385_S_CA)
     [  14.70431, 111.43862,  15.21419, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 385, "NCAC" ], // (385_S_N, 385_S_CA, 385_S_C)
     [  15.21419, 119.82201,  12.32404, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 385, "CACO" ], // (385_S_CA, 385_S_C, 385_S_O)
     [  15.21419, 117.59024,  13.37619, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 385, "CACN" ], // (385_S_CA, 385_S_C, 386_R_N)
     [  13.37619, 121.70888,  14.61551, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 385, "CNCA" ], // (385_S_C, 386_R_N, 386_R_CA)
     [  14.61551, 112.49691,  15.32733, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 386, "NCAC" ], // (386_R_N, 386_R_CA, 386_R_C)
     [  15.32733, 118.65334,  12.27969, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 386, "CACO" ], // (386_R_CA, 386_R_C, 386_R_O)
     [  15.32733, 119.47991,  13.36823, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 386, "CACN" ], // (386_R_CA, 386_R_C, 387_F_N)
     [  13.36823, 122.75572,  14.72974, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 386, "CNCA" ], // (386_R_C, 387_F_N, 387_F_CA)
     [  14.72974, 111.62817,  15.40729, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "F", 387, "NCAC" ], // (387_F_N, 387_F_CA, 387_F_C)
     [  15.40729, 119.33070,  12.34420, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "F", 387, "CACO" ], // (387_F_CA, 387_F_C, 387_F_O)
     [  15.40729, 118.20428,  13.45710, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "F", 387, "CACN" ], // (387_F_CA, 387_F_C, 388_R_N)
     [  13.45710, 123.89994,  14.79406, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "F", 387, "CNCA" ], // (387_F_C, 388_R_N, 388_R_CA)
     [  14.79406, 112.46577,  15.28015, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 388, "NCAC" ], // (388_R_N, 388_R_CA, 388_R_C)
     [  15.28015, 119.51730,  12.27083, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 388, "CACO" ], // (388_R_CA, 388_R_C, 388_R_O)
     [  15.28015, 118.02992,  13.33160, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 388, "CACN" ], // (388_R_CA, 388_R_C, 389_E_N)
     [  13.33160, 124.13177,  14.74607, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 388, "CNCA" ], // (388_R_C, 389_E_N, 389_E_CA)
     [  14.74607, 111.78457,  15.30433, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 389, "NCAC" ], // (389_E_N, 389_E_CA, 389_E_C)
     [  15.30433, 121.32093,  12.30098, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 389, "CACO" ], // (389_E_CA, 389_E_C, 389_E_O)
     [  15.30433, 116.13080,  13.34922, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 389, "CACN" ], // (389_E_CA, 389_E_C, 390_L_N)
     [  13.34922, 123.87817,  14.54576, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 389, "CNCA" ], // (389_E_C, 390_L_N, 390_L_CA)
     [  14.54576, 112.96798,  15.34093, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 390, "NCAC" ], // (390_L_N, 390_L_CA, 390_L_C)
     [  15.34093, 117.57291,  12.28582, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 390, "CACO" ], // (390_L_CA, 390_L_C, 390_L_O)
     [  15.34093, 120.55990,  13.39197, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 390, "CACN" ], // (390_L_CA, 390_L_C, 391_K_N)
     [  13.39197, 124.65044,  14.74654, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 390, "CNCA" ], // (390_L_C, 391_K_N, 391_K_CA)
     [  14.74654, 110.81253,  15.33904, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 391, "NCAC" ], // (391_K_N, 391_K_CA, 391_K_C)
     [  15.33904, 120.76490,  12.30558, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 391, "CACO" ], // (391_K_CA, 391_K_C, 391_K_O)
     [  15.33904, 116.29362,  13.36263, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 391, "CACN" ], // (391_K_CA, 391_K_C, 392_L_N)
     [  13.36263, 123.62444,  14.63711, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 391, "CNCA" ], // (391_K_C, 392_L_N, 392_L_CA)
     [  14.63711, 110.53640,  15.36216, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 392, "NCAC" ], // (392_L_N, 392_L_CA, 392_L_C)
     [  15.36216, 120.52713,  12.28203, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 392, "CACO" ], // (392_L_CA, 392_L_C, 392_L_O)
     [  15.36216, 117.56726,  13.32223, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 392, "CACN" ], // (392_L_CA, 392_L_C, 393_Q_N)
     [  13.32223, 122.23965,  14.59685, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 392, "CNCA" ], // (392_L_C, 393_Q_N, 393_Q_CA)
     [  14.59685, 110.48146,  15.28592, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 393, "NCAC" ], // (393_Q_N, 393_Q_CA, 393_Q_C)
     [  15.28592, 121.03923,  12.30286, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 393, "CACO" ], // (393_Q_CA, 393_Q_C, 393_Q_O)
     [  15.28592, 117.11184,  13.35620, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 393, "CACN" ], // (393_Q_CA, 393_Q_C, 394_H_N)
     [  13.35620, 122.28302,  14.66303, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 393, "CNCA" ], // (393_Q_C, 394_H_N, 394_H_CA)
     [  14.66303, 112.47921,  15.39696, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 394, "NCAC" ], // (394_H_N, 394_H_CA, 394_H_C)
     [  15.39696, 120.33001,  12.26045, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 394, "CACO" ], // (394_H_CA, 394_H_C, 394_H_O)
     [  15.39696, 118.15329,  13.40364, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 394, "CACN" ], // (394_H_CA, 394_H_C, 395_K_N)
     [  13.40364, 122.91256,  14.73148, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 394, "CNCA" ], // (394_H_C, 395_K_N, 395_K_CA)
     [  14.73148, 113.00955,  15.36996, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 395, "NCAC" ], // (395_K_N, 395_K_CA, 395_K_C)
     [  15.36996, 119.26910,  12.25107, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 395, "CACO" ], // (395_K_CA, 395_K_C, 395_K_O)
     [  15.36996, 118.57024,  13.33903, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 395, "CACN" ], // (395_K_CA, 395_K_C, 396_E_N)
     [  13.33903, 123.08374,  14.70632, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 395, "CNCA" ], // (395_K_C, 396_E_N, 396_E_CA)
     [  14.70632, 111.07781,  15.35754, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 396, "NCAC" ], // (396_E_N, 396_E_CA, 396_E_C)
     [  15.35754, 119.71492,  12.30114, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 396, "CACO" ], // (396_E_CA, 396_E_C, 396_E_O)
     [  15.35754, 117.80047,  13.39835, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 396, "CACN" ], // (396_E_CA, 396_E_C, 397_Y_N)
     [  13.39835, 124.12649,  14.67257, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 396, "CNCA" ], // (396_E_C, 397_Y_N, 397_Y_CA)
     [  14.67257, 111.77342,  15.36293, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Y", 397, "NCAC" ], // (397_Y_N, 397_Y_CA, 397_Y_C)
     [  15.36293, 119.28443,  12.31353, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Y", 397, "CACO" ], // (397_Y_CA, 397_Y_C, 397_Y_O)
     [  15.36293, 118.02129,  13.39701, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Y", 397, "CACN" ], // (397_Y_CA, 397_Y_C, 398_L_N)
     [  13.39701, 123.49570,  14.68995, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Y", 397, "CNCA" ], // (397_Y_C, 398_L_N, 398_L_CA)
     [  14.68995, 111.93172,  15.31590, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 398, "NCAC" ], // (398_L_N, 398_L_CA, 398_L_C)
     [  15.31590, 120.06951,  12.31515, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 398, "CACO" ], // (398_L_CA, 398_L_C, 398_L_O)
     [  15.31590, 117.63668,  13.35917, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 398, "CACN" ], // (398_L_CA, 398_L_C, 399_C_N)
     [  13.35917, 122.10328,  14.57964, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 398, "CNCA" ], // (398_L_C, 399_C_N, 399_C_CA)
     [  14.57964, 111.65043,  15.35652, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "C", 399, "NCAC" ], // (399_C_N, 399_C_CA, 399_C_C)
     [  15.35652, 118.72727,  12.25978, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "C", 399, "CACO" ], // (399_C_CA, 399_C_C, 399_C_O)
     [  15.35652, 119.90494,  13.39369, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "C", 399, "CACN" ], // (399_C_CA, 399_C_C, 400_V_N)
     [  13.39369, 122.46226,  14.80024, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "C", 399, "CNCA" ], // (399_C_C, 400_V_N, 400_V_CA)
     [  14.80024, 111.23604,  15.37995, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 400, "NCAC" ], // (400_V_N, 400_V_CA, 400_V_C)
     [  15.37995, 119.91946,  12.29889, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 400, "CACO" ], // (400_V_CA, 400_V_C, 400_V_O)
     [  15.37995, 118.56200,  13.39023, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 400, "CACN" ], // (400_V_CA, 400_V_C, 401_K_N)
     [  13.39023, 124.71490,  14.73366, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 400, "CNCA" ], // (400_V_C, 401_K_N, 401_K_CA)
     [  14.73366, 112.19233,  15.37116, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 401, "NCAC" ], // (401_K_N, 401_K_CA, 401_K_C)
     [  15.37116, 119.83024,  12.28528, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 401, "CACO" ], // (401_K_CA, 401_K_C, 401_K_O)
     [  15.37116, 118.41941,  13.33131, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 401, "CACN" ], // (401_K_CA, 401_K_C, 402_A_N)
     [  13.33131, 122.99704,  14.64880, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 401, "CNCA" ], // (401_K_C, 402_A_N, 402_A_CA)
     [  14.64880, 112.35370,  15.43705, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 402, "NCAC" ], // (402_A_N, 402_A_CA, 402_A_C)
     [  15.43705, 119.85398,  12.30443, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 402, "CACO" ], // (402_A_CA, 402_A_C, 402_A_O)
     [  15.43705, 118.75415,  13.40696, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 402, "CACN" ], // (402_A_CA, 402_A_C, 403_M_N)
     [  13.40696, 122.67567,  14.67739, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 402, "CNCA" ], // (402_A_C, 403_M_N, 403_M_CA)
     [  14.67739, 111.56266,  15.29674, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 403, "NCAC" ], // (403_M_N, 403_M_CA, 403_M_C)
     [  15.29674, 118.03315,  12.23071, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 403, "CACO" ], // (403_M_CA, 403_M_C, 403_M_O)
     [  15.29674, 120.04430,  13.39849, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 403, "CACN" ], // (403_M_CA, 403_M_C, 404_I_N)
     [  13.39849, 123.10792,  14.77237, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 403, "CNCA" ], // (403_M_C, 404_I_N, 404_I_CA)
     [  14.77237, 111.22722,  15.39404, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 404, "NCAC" ], // (404_I_N, 404_I_CA, 404_I_C)
     [  15.39404, 119.18086,  12.29600, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 404, "CACO" ], // (404_I_CA, 404_I_C, 404_I_O)
     [  15.39404, 118.39647,  13.40322, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 404, "CACN" ], // (404_I_CA, 404_I_C, 405_L_N)
     [  13.40322, 122.95212,  14.72074, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 404, "CNCA" ], // (404_I_C, 405_L_N, 405_L_CA)
     [  14.72074, 112.75272,  15.37739, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 405, "NCAC" ], // (405_L_N, 405_L_CA, 405_L_C)
     [  15.37739, 120.56420,  12.29315, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 405, "CACO" ], // (405_L_CA, 405_L_C, 405_L_O)
     [  15.37739, 117.32618,  13.34551, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 405, "CACN" ], // (405_L_CA, 405_L_C, 406_L_N)
     [  13.34551, 123.94597,  14.66903, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 405, "CNCA" ], // (405_L_C, 406_L_N, 406_L_CA)
     [  14.66903, 113.52252,  15.42482, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 406, "NCAC" ], // (406_L_N, 406_L_CA, 406_L_C)
     [  15.42482, 120.04709,  12.27996, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 406, "CACO" ], // (406_L_CA, 406_L_C, 406_L_O)
     [  15.42482, 118.75394,  13.34067, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 406, "CACN" ], // (406_L_CA, 406_L_C, 407_N_N)
     [  13.34067, 125.50072,  14.68041, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 406, "CNCA" ], // (406_L_C, 407_N_N, 407_N_CA)
     [  14.68041, 113.23031,  15.46034, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 407, "NCAC" ], // (407_N_N, 407_N_CA, 407_N_C)
     [  15.46034, 121.80139,  12.35021, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 407, "CACO" ], // (407_N_CA, 407_N_C, 407_N_O)
     [  15.46034, 117.67069,  13.36598, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 407, "CACN" ], // (407_N_CA, 407_N_C, 408_S_N)
     [  13.36598, 121.49464,  14.65896, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 407, "CNCA" ], // (407_N_C, 408_S_N, 408_S_CA)
     [  14.65896, 110.03622,  15.25247, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 408, "NCAC" ], // (408_S_N, 408_S_CA, 408_S_C)
     [  15.25247, 121.02747,  12.31549, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 408, "CACO" ], // (408_S_CA, 408_S_C, 408_S_O)
     [  15.25247, 117.20425,  13.31633, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 408, "CACN" ], // (408_S_CA, 408_S_C, 409_S_N)
     [  13.31633, 122.57328,  14.59762, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 408, "CNCA" ], // (408_S_C, 409_S_N, 409_S_CA)
     [  14.59762, 110.43103,  15.34231, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 409, "NCAC" ], // (409_S_N, 409_S_CA, 409_S_C)
     [  15.34231, 122.42999,  12.36635, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 409, "CACO" ], // (409_S_CA, 409_S_C, 409_S_O)
     [  15.34231, 116.43633,  13.34881, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 409, "CACN" ], // (409_S_CA, 409_S_C, 410_M_N)
     [  13.34881, 122.76287,  14.60522, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 409, "CNCA" ], // (409_S_C, 410_M_N, 410_M_CA)
     [  14.60522, 109.93078,  15.31161, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 410, "NCAC" ], // (410_M_N, 410_M_CA, 410_M_C)
     [  15.31161, 122.18223,  12.27597, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 410, "CACO" ], // (410_M_CA, 410_M_C, 410_M_O)
     [  15.31161, 115.89638,  13.33849, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 410, "CACN" ], // (410_M_CA, 410_M_C, 411_Y_N)
     [  13.33849, 121.99965,  14.58412, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 410, "CNCA" ], // (410_M_C, 411_Y_N, 411_Y_CA)
     [  14.58412, 108.06186,  15.39926, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Y", 411, "NCAC" ], // (411_Y_N, 411_Y_CA, 411_Y_C)
     [  15.39926, 117.78337,  12.39097, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Y", 411, "CACO" ], // (411_Y_CA, 411_Y_C, 411_Y_O)
     [  15.39926, 119.68932,  13.44064, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Y", 411, "CACN" ], // (411_Y_CA, 411_Y_C, 412_P_N)
     [  13.44064, 122.93676,  14.55644, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Y", 411, "CNCA" ], // (411_Y_C, 412_P_N, 412_P_CA)
     [  14.55644, 108.84612,  15.29873, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 412, "NCAC" ], // (412_P_N, 412_P_CA, 412_P_C)
     [  15.29873, 121.82685,  12.32537, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 412, "CACO" ], // (412_P_CA, 412_P_C, 412_P_O)
     [  15.29873, 115.44778,  13.29205, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 412, "CACN" ], // (412_P_CA, 412_P_C, 413_L_N)
     [  13.29205, 124.42973,  14.52405, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 412, "CNCA" ], // (412_P_C, 413_L_N, 413_L_CA)
     [  14.52405, 105.92054,  15.33399, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 413, "NCAC" ], // (413_L_N, 413_L_CA, 413_L_C)
     [  15.33399, 118.79475,  12.25581, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 413, "CACO" ], // (413_L_CA, 413_L_C, 413_L_O)
     [  15.33399, 116.05762,  13.42467, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 413, "CACN" ], // (413_L_CA, 413_L_C, 414_V_N)
     [  13.42467, 129.34848,  14.86316, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 413, "CNCA" ], // (413_L_C, 414_V_N, 414_V_CA)
     [  14.86316, 116.31994,  15.48178, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 414, "NCAC" ], // (414_V_N, 414_V_CA, 414_V_C)
     [  15.48178, 123.70504,  12.18537, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 414, "CACO" ], // (414_V_CA, 414_V_C, 414_V_O)
     [  15.48178, 114.79732,  13.41507, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 414, "CACN" ], // (414_V_CA, 414_V_C, 415_T_N)
     [  13.41507, 127.99433,  14.85407, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 414, "CNCA" ], // (414_V_C, 415_T_N, 415_T_CA)
     [  14.85407, 118.74063,  15.41634, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 415, "NCAC" ], // (415_T_N, 415_T_CA, 415_T_C)
     [  15.41634, 122.96038,  12.29090, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 415, "CACO" ], // (415_T_CA, 415_T_C, 415_T_O)
     [  15.41634, 114.39132,  13.39238, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 415, "CACN" ], // (415_T_CA, 415_T_C, 416_A_N)
     [  13.39238, 128.65559,  14.70920, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 415, "CNCA" ], // (415_T_C, 416_A_N, 416_A_CA)
     [  14.70920, 118.79895,  15.53836, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 416, "NCAC" ], // (416_A_N, 416_A_CA, 416_A_C)
     [  15.53836, 119.19596,  12.31492, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 416, "CACO" ], // (416_A_CA, 416_A_C, 416_A_O)
     [  15.53836, 119.89741,  13.35483, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 416, "CACN" ], // (416_A_CA, 416_A_C, 417_T_N)
     [  13.35483, 124.71987,  14.73968, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 416, "CNCA" ], // (416_A_C, 417_T_N, 417_T_CA)
     [  14.73968, 113.68882,  15.26785, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 417, "NCAC" ], // (417_T_N, 417_T_CA, 417_T_C)
     [  15.26785, 117.38139,  12.22941, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 417, "CACO" ], // (417_T_CA, 417_T_C, 417_T_O)
     [  15.26785, 120.04015,  13.40456, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 417, "CACN" ], // (417_T_CA, 417_T_C, 418_Q_N)
     [  13.40456, 122.35436,  14.76233, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 417, "CNCA" ], // (417_T_C, 418_Q_N, 418_Q_CA)
     [  14.76233, 111.43500,  15.37162, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 418, "NCAC" ], // (418_Q_N, 418_Q_CA, 418_Q_C)
     [  15.37162, 119.55272,  12.30746, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 418, "CACO" ], // (418_Q_CA, 418_Q_C, 418_Q_O)
     [  15.37162, 117.68310,  13.40698, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 418, "CACN" ], // (418_Q_CA, 418_Q_C, 419_D_N)
     [  13.40698, 124.57704,  14.74676, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 418, "CNCA" ], // (418_Q_C, 419_D_N, 419_D_CA)
     [  14.74676, 113.80365,  15.41052, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 419, "NCAC" ], // (419_D_N, 419_D_CA, 419_D_C)
     [  15.41052, 120.65084,  12.30293, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 419, "CACO" ], // (419_D_CA, 419_D_C, 419_D_O)
     [  15.41052, 118.16625,  13.38473, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 419, "CACN" ], // (419_D_CA, 419_D_C, 420_A_N)
     [  13.38473, 123.77543,  14.63216, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 419, "CNCA" ], // (419_D_C, 420_A_N, 420_A_CA)
     [  14.63216, 114.05267,  15.39818, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 420, "NCAC" ], // (420_A_N, 420_A_CA, 420_A_C)
     [  15.39818, 119.46696,  12.33279, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 420, "CACO" ], // (420_A_CA, 420_A_C, 420_A_O)
     [  15.39818, 118.69867,  13.33790, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 420, "CACN" ], // (420_A_CA, 420_A_C, 421_D_N)
     [  13.33790, 122.98426,  14.66597, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 420, "CNCA" ], // (420_A_C, 421_D_N, 421_D_CA)
     [  14.66597, 110.13058,  15.29748, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 421, "NCAC" ], // (421_D_N, 421_D_CA, 421_D_C)
     [  15.29748, 119.35670,  12.32371, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 421, "CACO" ], // (421_D_CA, 421_D_C, 421_D_O)
     [  15.29748, 117.95764,  13.37721, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 421, "CACN" ], // (421_D_CA, 421_D_C, 422_S_N)
     [  13.37721, 123.34597,  14.70115, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 421, "CNCA" ], // (421_D_C, 422_S_N, 422_S_CA)
     [  14.70115, 112.60091,  15.30065, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 422, "NCAC" ], // (422_S_N, 422_S_CA, 422_S_C)
     [  15.30065, 119.42786,  12.29234, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 422, "CACO" ], // (422_S_CA, 422_S_C, 422_S_O)
     [  15.30065, 118.78241,  13.40345, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 422, "CACN" ], // (422_S_CA, 422_S_C, 423_S_N)
     [  13.40345, 122.75786,  14.73320, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 422, "CNCA" ], // (422_S_C, 423_S_N, 423_S_CA)
     [  14.73320, 112.36906,  15.29337, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 423, "NCAC" ], // (423_S_N, 423_S_CA, 423_S_C)
     [  15.29337, 120.14112,  12.29783, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 423, "CACO" ], // (423_S_CA, 423_S_C, 423_S_O)
     [  15.29337, 118.05381,  13.43632, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 423, "CACN" ], // (423_S_CA, 423_S_C, 424_R_N)
     [  13.43632, 122.29052,  14.69732, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 423, "CNCA" ], // (423_S_C, 424_R_N, 424_R_CA)
     [  14.69732, 110.98668,  15.19918, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 424, "NCAC" ], // (424_R_N, 424_R_CA, 424_R_C)
     [  15.19918, 118.37775,  12.24542, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 424, "CACO" ], // (424_R_CA, 424_R_C, 424_R_O)
     [  15.19918, 118.45663,  13.38235, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 424, "CACN" ], // (424_R_CA, 424_R_C, 425_K_N)
     [  13.38235, 123.00514,  14.70805, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 424, "CNCA" ], // (424_R_C, 425_K_N, 425_K_CA)
     [  14.70805, 110.98465,  15.24200, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 425, "NCAC" ], // (425_K_N, 425_K_CA, 425_K_C)
     [  15.24200, 118.58575,  12.23070, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 425, "CACO" ], // (425_K_CA, 425_K_C, 425_K_O)
     [  15.24200, 118.35256,  13.33704, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 425, "CACN" ], // (425_K_CA, 425_K_C, 426_L_N)
     [  13.33704, 123.42592,  14.64533, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 425, "CNCA" ], // (425_K_C, 426_L_N, 426_L_CA)
     [  14.64533, 112.66118,  15.39479, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 426, "NCAC" ], // (426_L_N, 426_L_CA, 426_L_C)
     [  15.39479, 119.64834,  12.28440, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 426, "CACO" ], // (426_L_CA, 426_L_C, 426_L_O)
     [  15.39479, 118.30830,  13.37291, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 426, "CACN" ], // (426_L_CA, 426_L_C, 427_A_N)
     [  13.37291, 122.53198,  14.62993, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 426, "CNCA" ], // (426_L_C, 427_A_N, 427_A_CA)
     [  14.62993, 112.68932,  15.42166, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 427, "NCAC" ], // (427_A_N, 427_A_CA, 427_A_C)
     [  15.42166, 119.31256,  12.27505, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 427, "CACO" ], // (427_A_CA, 427_A_C, 427_A_O)
     [  15.42166, 119.10545,  13.38342, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 427, "CACN" ], // (427_A_CA, 427_A_C, 428_H_N)
     [  13.38342, 122.75673,  14.68015, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 427, "CNCA" ], // (427_A_C, 428_H_N, 428_H_CA)
     [  14.68015, 111.92216,  15.39462, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 428, "NCAC" ], // (428_H_N, 428_H_CA, 428_H_C)
     [  15.39462, 118.59446,  12.28970, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 428, "CACO" ], // (428_H_CA, 428_H_C, 428_H_O)
     [  15.39462, 119.07293,  13.42579, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 428, "CACN" ], // (428_H_CA, 428_H_C, 429_L_N)
     [  13.42579, 122.72058,  14.73274, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 428, "CNCA" ], // (428_H_C, 429_L_N, 429_L_CA)
     [  14.73274, 111.46095,  15.37289, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 429, "NCAC" ], // (429_L_N, 429_L_CA, 429_L_C)
     [  15.37289, 120.65833,  12.29729, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 429, "CACO" ], // (429_L_CA, 429_L_C, 429_L_O)
     [  15.37289, 117.33011,  13.37993, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 429, "CACN" ], // (429_L_CA, 429_L_C, 430_L_N)
     [  13.37993, 123.19751,  14.63166, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 429, "CNCA" ], // (429_L_C, 430_L_N, 430_L_CA)
     [  14.63166, 111.85031,  15.37430, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 430, "NCAC" ], // (430_L_N, 430_L_CA, 430_L_C)
     [  15.37430, 119.41359,  12.31831, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 430, "CACO" ], // (430_L_CA, 430_L_C, 430_L_O)
     [  15.37430, 118.25346,  13.37005, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 430, "CACN" ], // (430_L_CA, 430_L_C, 431_N_N)
     [  13.37005, 123.09319,  14.63401, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 430, "CNCA" ], // (430_L_C, 431_N_N, 431_N_CA)
     [  14.63401, 111.06159,  15.34697, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 431, "NCAC" ], // (431_N_N, 431_N_CA, 431_N_C)
     [  15.34697, 118.60494,  12.22866, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 431, "CACO" ], // (431_N_CA, 431_N_C, 431_N_O)
     [  15.34697, 119.27406,  13.36789, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 431, "CACN" ], // (431_N_CA, 431_N_C, 432_A_N)
     [  13.36789, 122.39020,  14.65637, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 431, "CNCA" ], // (431_N_C, 432_A_N, 432_A_CA)
     [  14.65637, 113.06581,  15.41661, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 432, "NCAC" ], // (432_A_N, 432_A_CA, 432_A_C)
     [  15.41661, 120.29354,  12.31793, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 432, "CACO" ], // (432_A_CA, 432_A_C, 432_A_O)
     [  15.41661, 118.54471,  13.39735, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 432, "CACN" ], // (432_A_CA, 432_A_C, 433_V_N)
     [  13.39735, 122.42603,  14.77103, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 432, "CNCA" ], // (432_A_C, 433_V_N, 433_V_CA)
     [  14.77103, 111.79570,  15.38535, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 433, "NCAC" ], // (433_V_N, 433_V_CA, 433_V_C)
     [  15.38535, 119.89885,  12.32046, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 433, "CACO" ], // (433_V_CA, 433_V_C, 433_V_O)
     [  15.38535, 118.57985,  13.34988, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 433, "CACN" ], // (433_V_CA, 433_V_C, 434_T_N)
     [  13.34988, 123.20744,  14.74832, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 433, "CNCA" ], // (433_V_C, 434_T_N, 434_T_CA)
     [  14.74832, 111.26559,  15.32308, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "T", 434, "NCAC" ], // (434_T_N, 434_T_CA, 434_T_C)
     [  15.32308, 119.53617,  12.29159, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "T", 434, "CACO" ], // (434_T_CA, 434_T_C, 434_T_O)
     [  15.32308, 118.51912,  13.33457, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "T", 434, "CACN" ], // (434_T_CA, 434_T_C, 435_D_N)
     [  13.33457, 124.52102,  14.65822, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "T", 434, "CNCA" ], // (434_T_C, 435_D_N, 435_D_CA)
     [  14.65822, 111.34349,  15.38761, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 435, "NCAC" ], // (435_D_N, 435_D_CA, 435_D_C)
     [  15.38761, 120.29415,  12.32303, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 435, "CACO" ], // (435_D_CA, 435_D_C, 435_D_O)
     [  15.38761, 118.02091,  13.38874, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 435, "CACN" ], // (435_D_CA, 435_D_C, 436_A_N)
     [  13.38874, 123.93309,  14.62333, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 435, "CNCA" ], // (435_D_C, 436_A_N, 436_A_CA)
     [  14.62333, 112.52921,  15.36397, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 436, "NCAC" ], // (436_A_N, 436_A_CA, 436_A_C)
     [  15.36397, 119.52465,  12.32830, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 436, "CACO" ], // (436_A_CA, 436_A_C, 436_A_O)
     [  15.36397, 118.68691,  13.36303, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 436, "CACN" ], // (436_A_CA, 436_A_C, 437_L_N)
     [  13.36303, 122.51667,  14.65026, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 436, "CNCA" ], // (436_A_C, 437_L_N, 437_L_CA)
     [  14.65026, 111.75629,  15.37942, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 437, "NCAC" ], // (437_L_N, 437_L_CA, 437_L_C)
     [  15.37942, 120.04162,  12.31085, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 437, "CACO" ], // (437_L_CA, 437_L_C, 437_L_O)
     [  15.37942, 118.17116,  13.35265, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 437, "CACN" ], // (437_L_CA, 437_L_C, 438_V_N)
     [  13.35265, 123.58738,  14.73977, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 437, "CNCA" ], // (437_L_C, 438_V_N, 438_V_CA)
     [  14.73977, 111.10500,  15.31940, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 438, "NCAC" ], // (438_V_N, 438_V_CA, 438_V_C)
     [  15.31940, 119.64446,  12.31277, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 438, "CACO" ], // (438_V_CA, 438_V_C, 438_V_O)
     [  15.31940, 118.42529,  13.35092, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 438, "CACN" ], // (438_V_CA, 438_V_C, 439_W_N)
     [  13.35092, 123.72502,  14.62908, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 438, "CNCA" ], // (438_V_C, 439_W_N, 439_W_CA)
     [  14.62908, 112.37905,  15.39995, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "W", 439, "NCAC" ], // (439_W_N, 439_W_CA, 439_W_C)
     [  15.39995, 119.23496,  12.32044, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "W", 439, "CACO" ], // (439_W_CA, 439_W_C, 439_W_O)
     [  15.39995, 119.13410,  13.41537, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "W", 439, "CACN" ], // (439_W_CA, 439_W_C, 440_V_N)
     [  13.41537, 123.41143,  14.80896, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "W", 439, "CNCA" ], // (439_W_C, 440_V_N, 440_V_CA)
     [  14.80896, 111.04826,  15.33239, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 440, "NCAC" ], // (440_V_N, 440_V_CA, 440_V_C)
     [  15.33239, 120.58100,  12.34103, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 440, "CACO" ], // (440_V_CA, 440_V_C, 440_V_O)
     [  15.33239, 117.40688,  13.34947, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 440, "CACN" ], // (440_V_CA, 440_V_C, 441_I_N)
     [  13.34947, 123.65283,  14.68167, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 440, "CNCA" ], // (440_V_C, 441_I_N, 441_I_CA)
     [  14.68167, 110.57656,  15.38325, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 441, "NCAC" ], // (441_I_N, 441_I_CA, 441_I_C)
     [  15.38325, 119.99647,  12.32471, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 441, "CACO" ], // (441_I_CA, 441_I_C, 441_I_O)
     [  15.38325, 118.58046,  13.36053, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 441, "CACN" ], // (441_I_CA, 441_I_C, 442_A_N)
     [  13.36053, 123.41467,  14.62806, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 441, "CNCA" ], // (441_I_C, 442_A_N, 442_A_CA)
     [  14.62806, 112.60018,  15.42159, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 442, "NCAC" ], // (442_A_N, 442_A_CA, 442_A_C)
     [  15.42159, 120.36650,  12.33021, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 442, "CACO" ], // (442_A_CA, 442_A_C, 442_A_O)
     [  15.42159, 118.22929,  13.41699, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 442, "CACN" ], // (442_A_CA, 442_A_C, 443_K_N)
     [  13.41699, 124.02670,  14.74688, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 442, "CNCA" ], // (442_A_C, 443_K_N, 443_K_CA)
     [  14.74688, 114.10542,  15.33067, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 443, "NCAC" ], // (443_K_N, 443_K_CA, 443_K_C)
     [  15.33067, 118.45141,  12.27179, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 443, "CACO" ], // (443_K_CA, 443_K_C, 443_K_O)
     [  15.33067, 119.29292,  13.35044, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 443, "CACN" ], // (443_K_CA, 443_K_C, 444_S_N)
     [  13.35044, 123.43193,  14.73156, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 443, "CNCA" ], // (443_K_C, 444_S_N, 444_S_CA)
     [  14.73156, 114.24368,  15.39390, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 444, "NCAC" ], // (444_S_N, 444_S_CA, 444_S_C)
     [  15.39390, 119.44272,  12.27215, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 444, "CACO" ], // (444_S_CA, 444_S_C, 444_S_O)
     [  15.39390, 119.26455,  13.40904, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 444, "CACN" ], // (444_S_CA, 444_S_C, 445_G_N)
     [  13.40904, 123.64116,  14.72750, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 444, "CNCA" ], // (444_S_C, 445_G_N, 445_G_CA)
     [  14.72750, 115.71146,  15.35491, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 445, "NCAC" ], // (445_G_N, 445_G_CA, 445_G_C)
     [  15.35491, 119.93964,  12.28639, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 445, "CACO" ], // (445_G_CA, 445_G_C, 445_G_O)
     [  15.35491, 118.31017,  13.35790, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 445, "CACN" ], // (445_G_CA, 445_G_C, 446_I_N)
     [  13.35790, 124.00850,  14.70110, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 445, "CNCA" ], // (445_G_C, 446_I_N, 446_I_CA)
     [  14.70110, 111.81581,  15.39826, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "I", 446, "NCAC" ], // (446_I_N, 446_I_CA, 446_I_C)
     [  15.39826, 121.14931,  12.34178, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "I", 446, "CACO" ], // (446_I_CA, 446_I_C, 446_I_O)
     [  15.39826, 117.22980,  13.36283, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "I", 446, "CACN" ], // (446_I_CA, 446_I_C, 447_S_N)
     [  13.36283, 123.14938,  14.74047, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "I", 446, "CNCA" ], // (446_I_C, 447_S_N, 447_S_CA)
     [  14.74047, 112.64209,  15.34478, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 447, "NCAC" ], // (447_S_N, 447_S_CA, 447_S_C)
     [  15.34478, 120.25313,  12.27285, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 447, "CACO" ], // (447_S_CA, 447_S_C, 447_S_O)
     [  15.34478, 118.38356,  13.38472, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 447, "CACN" ], // (447_S_CA, 447_S_C, 448_S_N)
     [  13.38472, 123.57971,  14.78034, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 447, "CNCA" ], // (447_S_C, 448_S_N, 448_S_CA)
     [  14.78034, 113.75369,  15.42523, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 448, "NCAC" ], // (448_S_N, 448_S_CA, 448_S_C)
     [  15.42523, 120.57266,  12.30793, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 448, "CACO" ], // (448_S_CA, 448_S_C, 448_S_O)
     [  15.42523, 118.29892,  13.39687, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 448, "CACN" ], // (448_S_CA, 448_S_C, 449_Q_N)
     [  13.39687, 122.81379,  14.70057, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 448, "CNCA" ], // (448_S_C, 449_Q_N, 449_Q_CA)
     [  14.70057, 112.34812,  15.40413, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 449, "NCAC" ], // (449_Q_N, 449_Q_CA, 449_Q_C)
     [  15.40413, 120.42356,  12.33790, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 449, "CACO" ], // (449_Q_CA, 449_Q_C, 449_Q_O)
     [  15.40413, 117.73100,  13.37156, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 449, "CACN" ], // (449_Q_CA, 449_Q_C, 450_Q_N)
     [  13.37156, 123.36479,  14.63891, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 449, "CNCA" ], // (449_Q_C, 450_Q_N, 450_Q_CA)
     [  14.63891, 111.98596,  15.36341, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 450, "NCAC" ], // (450_Q_N, 450_Q_CA, 450_Q_C)
     [  15.36341, 118.25436,  12.25432, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 450, "CACO" ], // (450_Q_CA, 450_Q_C, 450_Q_O)
     [  15.36341, 119.90531,  13.38875, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 450, "CACN" ], // (450_Q_CA, 450_Q_C, 451_Q_N)
     [  13.38875, 123.19251,  14.71318, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 450, "CNCA" ], // (450_Q_C, 451_Q_N, 451_Q_CA)
     [  14.71318, 110.98452,  15.35243, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Q", 451, "NCAC" ], // (451_Q_N, 451_Q_CA, 451_Q_C)
     [  15.35243, 119.79480,  12.30349, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Q", 451, "CACO" ], // (451_Q_CA, 451_Q_C, 451_Q_O)
     [  15.35243, 118.49192,  13.39202, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Q", 451, "CACN" ], // (451_Q_CA, 451_Q_C, 452_S_N)
     [  13.39202, 122.42990,  14.72482, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Q", 451, "CNCA" ], // (451_Q_C, 452_S_N, 452_S_CA)
     [  14.72482, 112.33906,  15.33431, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 452, "NCAC" ], // (452_S_N, 452_S_CA, 452_S_C)
     [  15.33431, 120.15738,  12.30584, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 452, "CACO" ], // (452_S_CA, 452_S_C, 452_S_O)
     [  15.33431, 118.06486,  13.39100, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 452, "CACN" ], // (452_S_CA, 452_S_C, 453_M_N)
     [  13.39100, 122.86349,  14.69442, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 452, "CNCA" ], // (452_S_C, 453_M_N, 453_M_CA)
     [  14.69442, 112.17607,  15.34361, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 453, "NCAC" ], // (453_M_N, 453_M_CA, 453_M_C)
     [  15.34361, 119.72933,  12.32035, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 453, "CACO" ], // (453_M_CA, 453_M_C, 453_M_O)
     [  15.34361, 118.37283,  13.38784, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 453, "CACN" ], // (453_M_CA, 453_M_C, 454_R_N)
     [  13.38784, 124.06564,  14.65642, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 453, "CNCA" ], // (453_M_C, 454_R_N, 454_R_CA)
     [  14.65642, 110.57383,  15.21527, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 454, "NCAC" ], // (454_R_N, 454_R_CA, 454_R_C)
     [  15.21527, 118.62708,  12.26530, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 454, "CACO" ], // (454_R_CA, 454_R_C, 454_R_O)
     [  15.21527, 117.50676,  13.32928, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 454, "CACN" ], // (454_R_CA, 454_R_C, 455_L_N)
     [  13.32928, 123.66413,  14.60754, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 454, "CNCA" ], // (454_R_C, 455_L_N, 455_L_CA)
     [  14.60754, 111.04986,  15.34589, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 455, "NCAC" ], // (455_L_N, 455_L_CA, 455_L_C)
     [  15.34589, 119.57679,  12.26217, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 455, "CACO" ], // (455_L_CA, 455_L_C, 455_L_O)
     [  15.34589, 118.31593,  13.36915, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 455, "CACN" ], // (455_L_CA, 455_L_C, 456_A_N)
     [  13.36915, 122.57966,  14.63038, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 455, "CNCA" ], // (455_L_C, 456_A_N, 456_A_CA)
     [  14.63038, 112.60459,  15.41939, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 456, "NCAC" ], // (456_A_N, 456_A_CA, 456_A_C)
     [  15.41939, 119.76556,  12.28708, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 456, "CACO" ], // (456_A_CA, 456_A_C, 456_A_O)
     [  15.41939, 118.59726,  13.39365, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 456, "CACN" ], // (456_A_CA, 456_A_C, 457_N_N)
     [  13.39365, 122.57359,  14.68370, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 456, "CNCA" ], // (456_A_C, 457_N_N, 457_N_CA)
     [  14.68370, 112.87858,  15.37170, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 457, "NCAC" ], // (457_N_N, 457_N_CA, 457_N_C)
     [  15.37170, 119.05253,  12.27509, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 457, "CACO" ], // (457_N_CA, 457_N_C, 457_N_O)
     [  15.37170, 118.46338,  13.36968, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 457, "CACN" ], // (457_N_CA, 457_N_C, 458_L_N)
     [  13.36968, 123.13330,  14.64911, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 457, "CNCA" ], // (457_N_C, 458_L_N, 458_L_CA)
     [  14.64911, 112.21739,  15.32856, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 458, "NCAC" ], // (458_L_N, 458_L_CA, 458_L_C)
     [  15.32856, 119.64576,  12.26469, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 458, "CACO" ], // (458_L_CA, 458_L_C, 458_L_O)
     [  15.32856, 118.40675,  13.36152, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 458, "CACN" ], // (458_L_CA, 458_L_C, 459_L_N)
     [  13.36152, 122.70391,  14.65780, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 458, "CNCA" ], // (458_L_C, 459_L_N, 459_L_CA)
     [  14.65780, 112.28308,  15.38961, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 459, "NCAC" ], // (459_L_N, 459_L_CA, 459_L_C)
     [  15.38961, 119.38328,  12.27776, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 459, "CACO" ], // (459_L_CA, 459_L_C, 459_L_O)
     [  15.38961, 118.98478,  13.38213, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 459, "CACN" ], // (459_L_CA, 459_L_C, 460_M_N)
     [  13.38213, 122.81526,  14.72230, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 459, "CNCA" ], // (459_L_C, 460_M_N, 460_M_CA)
     [  14.72230, 113.51951,  15.40058, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 460, "NCAC" ], // (460_M_N, 460_M_CA, 460_M_C)
     [  15.40058, 119.68766,  12.31477, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 460, "CACO" ], // (460_M_CA, 460_M_C, 460_M_O)
     [  15.40058, 119.53310,  13.39877, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 460, "CACN" ], // (460_M_CA, 460_M_C, 461_L_N)
     [  13.39877, 121.62912,  14.68495, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 460, "CNCA" ], // (460_M_C, 461_L_N, 461_L_CA)
     [  14.68495, 111.37063,  15.34764, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 461, "NCAC" ], // (461_L_N, 461_L_CA, 461_L_C)
     [  15.34764, 119.07693,  12.28546, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 461, "CACO" ], // (461_L_CA, 461_L_C, 461_L_O)
     [  15.34764, 118.36845,  13.40187, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 461, "CACN" ], // (461_L_CA, 461_L_C, 462_L_N)
     [  13.40187, 123.76537,  14.67660, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 461, "CNCA" ], // (461_L_C, 462_L_N, 462_L_CA)
     [  14.67660, 111.80430,  15.38167, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 462, "NCAC" ], // (462_L_N, 462_L_CA, 462_L_C)
     [  15.38167, 119.65441,  12.28474, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 462, "CACO" ], // (462_L_CA, 462_L_C, 462_L_O)
     [  15.38167, 119.01514,  13.41370, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 462, "CACN" ], // (462_L_CA, 462_L_C, 463_S_N)
     [  13.41370, 121.65130,  14.70525, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 462, "CNCA" ], // (462_L_C, 463_S_N, 463_S_CA)
     [  14.70525, 113.15188,  15.29913, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 463, "NCAC" ], // (463_S_N, 463_S_CA, 463_S_C)
     [  15.29913, 119.36116,  12.29666, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 463, "CACO" ], // (463_S_CA, 463_S_C, 463_S_O)
     [  15.29913, 118.39450,  13.37720, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 463, "CACN" ], // (463_S_CA, 463_S_C, 464_H_N)
     [  13.37720, 122.33048,  14.67382, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 463, "CNCA" ], // (463_S_C, 464_H_N, 464_H_CA)
     [  14.67382, 111.33682,  15.34976, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 464, "NCAC" ], // (464_H_N, 464_H_CA, 464_H_C)
     [  15.34976, 119.58861,  12.30449, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 464, "CACO" ], // (464_H_CA, 464_H_C, 464_H_O)
     [  15.34976, 117.57841,  13.40479, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 464, "CACN" ], // (464_H_CA, 464_H_C, 465_V_N)
     [  13.40479, 123.85412,  14.71778, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 464, "CNCA" ], // (464_H_C, 465_V_N, 465_V_CA)
     [  14.71778, 110.45376,  15.32556, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 465, "NCAC" ], // (465_V_N, 465_V_CA, 465_V_C)
     [  15.32556, 119.34766,  12.28670, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 465, "CACO" ], // (465_V_CA, 465_V_C, 465_V_O)
     [  15.32556, 118.62623,  13.42659, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 465, "CACN" ], // (465_V_CA, 465_V_C, 466_R_N)
     [  13.42659, 123.66558,  14.74182, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 465, "CNCA" ], // (465_V_C, 466_R_N, 466_R_CA)
     [  14.74182, 112.30709,  15.27727, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 466, "NCAC" ], // (466_R_N, 466_R_CA, 466_R_C)
     [  15.27727, 119.82406,  12.29187, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 466, "CACO" ], // (466_R_CA, 466_R_C, 466_R_O)
     [  15.27727, 117.96040,  13.32473, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 466, "CACN" ], // (466_R_CA, 466_R_C, 467_H_N)
     [  13.32473, 122.29262,  14.61701, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 466, "CNCA" ], // (466_R_C, 467_H_N, 467_H_CA)
     [  14.61701, 112.08227,  15.36591, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 467, "NCAC" ], // (467_H_N, 467_H_CA, 467_H_C)
     [  15.36591, 118.96999,  12.33270, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 467, "CACO" ], // (467_H_CA, 467_H_C, 467_H_O)
     [  15.36591, 118.10892,  13.38316, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 467, "CACN" ], // (467_H_CA, 467_H_C, 468_A_N)
     [  13.38316, 123.35263,  14.64512, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 467, "CNCA" ], // (467_H_C, 468_A_N, 468_A_CA)
     [  14.64512, 111.57639,  15.37697, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 468, "NCAC" ], // (468_A_N, 468_A_CA, 468_A_C)
     [  15.37697, 120.20205,  12.29215, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 468, "CACO" ], // (468_A_CA, 468_A_C, 468_A_O)
     [  15.37697, 117.52407,  13.36393, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 468, "CACN" ], // (468_A_CA, 468_A_C, 469_S_N)
     [  13.36393, 124.25541,  14.64338, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 468, "CNCA" ], // (468_A_C, 469_S_N, 469_S_CA)
     [  14.64338, 112.51250,  15.27370, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "S", 469, "NCAC" ], // (469_S_N, 469_S_CA, 469_S_C)
     [  15.27370, 119.34009,  12.27234, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "S", 469, "CACO" ], // (469_S_CA, 469_S_C, 469_S_O)
     [  15.27370, 118.76146,  13.37259, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "S", 469, "CACN" ], // (469_S_CA, 469_S_C, 470_N_N)
     [  13.37259, 122.15855,  14.65996, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "S", 469, "CNCA" ], // (469_S_C, 470_N_N, 470_N_CA)
     [  14.65996, 111.27896,  15.39470, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 470, "NCAC" ], // (470_N_N, 470_N_CA, 470_N_C)
     [  15.39470, 119.54666,  12.27454, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 470, "CACO" ], // (470_N_CA, 470_N_C, 470_N_O)
     [  15.39470, 117.97595,  13.40470, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 470, "CACN" ], // (470_N_CA, 470_N_C, 471_K_N)
     [  13.40470, 123.36864,  14.74140, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 470, "CNCA" ], // (470_N_C, 471_K_N, 471_K_CA)
     [  14.74140, 112.50884,  15.34996, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 471, "NCAC" ], // (471_K_N, 471_K_CA, 471_K_C)
     [  15.34996, 120.08450,  12.27990, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 471, "CACO" ], // (471_K_CA, 471_K_C, 471_K_O)
     [  15.34996, 117.95272,  13.35285, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 471, "CACN" ], // (471_K_CA, 471_K_C, 472_G_N)
     [  13.35285, 123.61797,  14.65997, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 471, "CNCA" ], // (471_K_C, 472_G_N, 472_G_CA)
     [  14.65997, 113.20710,  15.31122, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 472, "NCAC" ], // (472_G_N, 472_G_CA, 472_G_C)
     [  15.31122, 119.96334,  12.30000, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 472, "CACO" ], // (472_G_CA, 472_G_C, 472_G_O)
     [  15.31122, 117.53349,  13.36572, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "G", 472, "CACN" ], // (472_G_CA, 472_G_C, 473_M_N)
     [  13.36572, 124.22841,  14.65154, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "G", 472, "CNCA" ], // (472_G_C, 473_M_N, 473_M_CA)
     [  14.65154, 111.55179,  15.32142, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 473, "NCAC" ], // (473_M_N, 473_M_CA, 473_M_C)
     [  15.32142, 120.02494,  12.28172, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 473, "CACO" ], // (473_M_CA, 473_M_C, 473_M_O)
     [  15.32142, 117.40363,  13.36777, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 473, "CACN" ], // (473_M_CA, 473_M_C, 474_E_N)
     [  13.36777, 124.73806,  14.68492, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 473, "CNCA" ], // (473_M_C, 474_E_N, 474_E_CA)
     [  14.68492, 111.36005,  15.34352, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 474, "NCAC" ], // (474_E_N, 474_E_CA, 474_E_C)
     [  15.34352, 120.48246,  12.31402, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 474, "CACO" ], // (474_E_CA, 474_E_C, 474_E_O)
     [  15.34352, 117.83633,  13.35066, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 474, "CACN" ], // (474_E_CA, 474_E_C, 475_H_N)
     [  13.35066, 122.85218,  14.60007, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 474, "CNCA" ], // (474_E_C, 475_H_N, 475_H_CA)
     [  14.60007, 112.29079,  15.39202, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 475, "NCAC" ], // (475_H_N, 475_H_CA, 475_H_C)
     [  15.39202, 119.55431,  12.31229, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 475, "CACO" ], // (475_H_CA, 475_H_C, 475_H_O)
     [  15.39202, 118.95511,  13.38587, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 475, "CACN" ], // (475_H_CA, 475_H_C, 476_L_N)
     [  13.38587, 122.30376,  14.66886, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 475, "CNCA" ], // (475_H_C, 476_L_N, 476_L_CA)
     [  14.66886, 112.27007,  15.41109, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 476, "NCAC" ], // (476_L_N, 476_L_CA, 476_L_C)
     [  15.41109, 119.62312,  12.32676, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 476, "CACO" ], // (476_L_CA, 476_L_C, 476_L_O)
     [  15.41109, 118.38884,  13.38425, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 476, "CACN" ], // (476_L_CA, 476_L_C, 477_L_N)
     [  13.38425, 123.52499,  14.67915, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 476, "CNCA" ], // (476_L_C, 477_L_N, 477_L_CA)
     [  14.67915, 111.43235,  15.32532, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 477, "NCAC" ], // (477_L_N, 477_L_CA, 477_L_C)
     [  15.32532, 119.45389,  12.30904, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 477, "CACO" ], // (477_L_CA, 477_L_C, 477_L_O)
     [  15.32532, 118.15130,  13.33808, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 477, "CACN" ], // (477_L_CA, 477_L_C, 478_N_N)
     [  13.33808, 122.95137,  14.62481, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 477, "CNCA" ], // (477_L_C, 478_N_N, 478_N_CA)
     [  14.62481, 112.01413,  15.38566, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 478, "NCAC" ], // (478_N_N, 478_N_CA, 478_N_C)
     [  15.38566, 119.74628,  12.30753, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 478, "CACO" ], // (478_N_CA, 478_N_C, 478_N_O)
     [  15.38566, 118.57087,  13.36695, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 478, "CACN" ], // (478_N_CA, 478_N_C, 479_M_N)
     [  13.36695, 123.23289,  14.72153, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 478, "CNCA" ], // (478_N_C, 479_M_N, 479_M_CA)
     [  14.72153, 112.66139,  15.34928, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 479, "NCAC" ], // (479_M_N, 479_M_CA, 479_M_C)
     [  15.34928, 119.85954,  12.33589, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 479, "CACO" ], // (479_M_CA, 479_M_C, 479_M_O)
     [  15.34928, 117.71862,  13.41418, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 479, "CACN" ], // (479_M_CA, 479_M_C, 480_K_N)
     [  13.41418, 123.68865,  14.69657, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 479, "CNCA" ], // (479_M_C, 480_K_N, 480_K_CA)
     [  14.69657, 111.17320,  15.33527, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 480, "NCAC" ], // (480_K_N, 480_K_CA, 480_K_C)
     [  15.33527, 119.13974,  12.24113, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 480, "CACO" ], // (480_K_CA, 480_K_C, 480_K_O)
     [  15.33527, 119.07179,  13.34662, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 480, "CACN" ], // (480_K_CA, 480_K_C, 481_C_N)
     [  13.34662, 123.24804,  14.69708, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 480, "CNCA" ], // (480_K_C, 481_C_N, 481_C_CA)
     [  14.69708, 113.09594,  15.38069, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "C", 481, "NCAC" ], // (481_C_N, 481_C_CA, 481_C_C)
     [  15.38069, 119.71840,  12.29036, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "C", 481, "CACO" ], // (481_C_CA, 481_C_C, 481_C_O)
     [  15.38069, 118.05141,  13.38465, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "C", 481, "CACN" ], // (481_C_CA, 481_C_C, 482_K_N)
     [  13.38465, 124.41389,  14.71070, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "C", 481, "CNCA" ], // (481_C_C, 482_K_N, 482_K_CA)
     [  14.71070, 113.95499,  15.46240, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "K", 482, "NCAC" ], // (482_K_N, 482_K_CA, 482_K_C)
     [  15.46240, 118.72988,  12.25406, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "K", 482, "CACO" ], // (482_K_CA, 482_K_C, 482_K_O)
     [  15.46240, 119.60068,  13.39503, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "K", 482, "CACN" ], // (482_K_CA, 482_K_C, 483_N_N)
     [  13.39503, 126.15678,  14.73692, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "K", 482, "CNCA" ], // (482_K_C, 483_N_N, 483_N_CA)
     [  14.73692, 114.75125,  15.46373, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 483, "NCAC" ], // (483_N_N, 483_N_CA, 483_N_C)
     [  15.46373, 120.25987,  12.28948, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 483, "CACO" ], // (483_N_CA, 483_N_C, 483_N_O)
     [  15.46373, 118.57377,  13.40646, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 483, "CACN" ], // (483_N_CA, 483_N_C, 484_V_N)
     [  13.40646, 122.67714,  14.81456, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 483, "CNCA" ], // (483_N_C, 484_V_N, 484_V_CA)
     [  14.81456, 111.92490,  15.32208, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 484, "NCAC" ], // (484_V_N, 484_V_CA, 484_V_C)
     [  15.32208, 120.36733,  12.29688, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 484, "CACO" ], // (484_V_CA, 484_V_C, 484_V_O)
     [  15.32208, 117.48390,  13.37822, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 484, "CACN" ], // (484_V_CA, 484_V_C, 485_V_N)
     [  13.37822, 124.55297,  14.70106, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 484, "CNCA" ], // (484_V_C, 485_V_N, 485_V_CA)
     [  14.70106, 111.07611,  15.32234, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 485, "NCAC" ], // (485_V_N, 485_V_CA, 485_V_C)
     [  15.32234, 118.97158,  12.34770, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 485, "CACO" ], // (485_V_CA, 485_V_C, 485_V_O)
     [  15.32234, 119.09426,  13.50062, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 485, "CACN" ], // (485_V_CA, 485_V_C, 486_P_N)
     [  13.50062, 122.45622,  14.57605, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 485, "CNCA" ], // (485_V_C, 486_P_N, 486_P_CA)
     [  14.57605, 111.11109,  15.34227, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "P", 486, "NCAC" ], // (486_P_N, 486_P_CA, 486_P_C)
     [  15.34227, 121.72721,  12.33577, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "P", 486, "CACO" ], // (486_P_CA, 486_P_C, 486_P_O)
     [  15.34227, 115.87775,  13.30258, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "P", 486, "CACN" ], // (486_P_CA, 486_P_C, 487_V_N)
     [  13.30258, 123.89337,  14.65511, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "P", 486, "CNCA" ], // (486_P_C, 487_V_N, 487_V_CA)
     [  14.65511, 107.56199,  15.27611, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 487, "NCAC" ], // (487_V_N, 487_V_CA, 487_V_C)
     [  15.27611, 119.38298,  12.27142, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 487, "CACO" ], // (487_V_CA, 487_V_C, 487_V_O)
     [  15.27611, 117.05108,  13.33937, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 487, "CACN" ], // (487_V_CA, 487_V_C, 488_Y_N)
     [  13.33937, 124.53443,  14.68118, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 487, "CNCA" ], // (487_V_C, 488_Y_N, 488_Y_CA)
     [  14.68118, 111.69272,  15.35249, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "Y", 488, "NCAC" ], // (488_Y_N, 488_Y_CA, 488_Y_C)
     [  15.35249, 119.22452,  12.30153, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "Y", 488, "CACO" ], // (488_Y_CA, 488_Y_C, 488_Y_O)
     [  15.35249, 118.03395,  13.36395, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "Y", 488, "CACN" ], // (488_Y_CA, 488_Y_C, 489_D_N)
     [  13.36395, 123.89227,  14.73219, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "Y", 488, "CNCA" ], // (488_Y_C, 489_D_N, 489_D_CA)
     [  14.73219, 113.40413,  15.39284, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "D", 489, "NCAC" ], // (489_D_N, 489_D_CA, 489_D_C)
     [  15.39284, 119.98971,  12.30408, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "D", 489, "CACO" ], // (489_D_CA, 489_D_C, 489_D_O)
     [  15.39284, 118.39908,  13.39815, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "D", 489, "CACN" ], // (489_D_CA, 489_D_C, 490_L_N)
     [  13.39815, 122.80648,  14.68033, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "D", 489, "CNCA" ], // (489_D_C, 490_L_N, 490_L_CA)
     [  14.68033, 111.83243,  15.35935, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 490, "NCAC" ], // (490_L_N, 490_L_CA, 490_L_C)
     [  15.35935, 120.68966,  12.31341, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 490, "CACO" ], // (490_L_CA, 490_L_C, 490_L_O)
     [  15.35935, 117.98951,  13.36033, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 490, "CACN" ], // (490_L_CA, 490_L_C, 491_L_N)
     [  13.36033, 122.54170,  14.65464, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 490, "CNCA" ], // (490_L_C, 491_L_N, 491_L_CA)
     [  14.65464, 111.52544,  15.41483, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 491, "NCAC" ], // (491_L_N, 491_L_CA, 491_L_C)
     [  15.41483, 120.10659,  12.29997, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 491, "CACO" ], // (491_L_CA, 491_L_C, 491_L_O)
     [  15.41483, 118.00303,  13.34345, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 491, "CACN" ], // (491_L_CA, 491_L_C, 492_L_N)
     [  13.34345, 125.43790,  14.59921, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 491, "CNCA" ], // (491_L_C, 492_L_N, 492_L_CA)
     [  14.59921, 111.82863,  15.36463, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 492, "NCAC" ], // (492_L_N, 492_L_CA, 492_L_C)
     [  15.36463, 118.99588,  12.31903, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 492, "CACO" ], // (492_L_CA, 492_L_C, 492_L_O)
     [  15.36463, 118.64150,  13.34044, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 492, "CACN" ], // (492_L_CA, 492_L_C, 493_E_N)
     [  13.34044, 124.38252,  14.70125, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 492, "CNCA" ], // (492_L_C, 493_E_N, 493_E_CA)
     [  14.70125, 111.19744,  15.36383, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "E", 493, "NCAC" ], // (493_E_N, 493_E_CA, 493_E_C)
     [  15.36383, 120.49214,  12.36749, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "E", 493, "CACO" ], // (493_E_CA, 493_E_C, 493_E_O)
     [  15.36383, 117.03615,  13.36966, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "E", 493, "CACN" ], // (493_E_CA, 493_E_C, 494_M_N)
     [  13.36966, 124.35534,  14.69404, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "E", 493, "CNCA" ], // (493_E_C, 494_M_N, 494_M_CA)
     [  14.69404, 112.44978,  15.39581, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "M", 494, "NCAC" ], // (494_M_N, 494_M_CA, 494_M_C)
     [  15.39581, 120.79404,  12.32823, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "M", 494, "CACO" ], // (494_M_CA, 494_M_C, 494_M_O)
     [  15.39581, 117.75172,  13.35273, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "M", 494, "CACN" ], // (494_M_CA, 494_M_C, 495_L_N)
     [  13.35273, 124.97288,  14.64587, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "M", 494, "CNCA" ], // (494_M_C, 495_L_N, 495_L_CA)
     [  14.64587, 112.17509,  15.43029, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 495, "NCAC" ], // (495_L_N, 495_L_CA, 495_L_C)
     [  15.43029, 120.67181,  12.34021, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 495, "CACO" ], // (495_L_CA, 495_L_C, 495_L_O)
     [  15.43029, 118.09750,  13.37081, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 495, "CACN" ], // (495_L_CA, 495_L_C, 496_N_N)
     [  13.37081, 123.87945,  14.67641, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 495, "CNCA" ], // (495_L_C, 496_N_N, 496_N_CA)
     [  14.67641, 112.43648,  15.40751, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "N", 496, "NCAC" ], // (496_N_N, 496_N_CA, 496_N_C)
     [  15.40751, 119.94789,  12.36973, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "N", 496, "CACO" ], // (496_N_CA, 496_N_C, 496_N_O)
     [  15.40751, 118.02245,  13.37725, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "N", 496, "CACN" ], // (496_N_CA, 496_N_C, 497_A_N)
     [  13.37725, 123.96098,  14.62780, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "N", 496, "CNCA" ], // (496_N_C, 497_A_N, 497_A_CA)
     [  14.62780, 113.19279,  15.43550, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "A", 497, "NCAC" ], // (497_A_N, 497_A_CA, 497_A_C)
     [  15.43550, 120.06007,  12.34317, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "A", 497, "CACO" ], // (497_A_CA, 497_A_C, 497_A_O)
     [  15.43550, 118.63373,  13.39365, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "A", 497, "CACN" ], // (497_A_CA, 497_A_C, 498_H_N)
     [  13.39365, 122.89005,  14.67731, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "A", 497, "CNCA" ], // (497_A_C, 498_H_N, 498_H_CA)
     [  14.67731, 112.89204,  15.36526, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "H", 498, "NCAC" ], // (498_H_N, 498_H_CA, 498_H_C)
     [  15.36526, 119.30223,  12.33558, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "H", 498, "CACO" ], // (498_H_CA, 498_H_C, 498_H_O)
     [  15.36526, 119.02347,  13.39499, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "H", 498, "CACN" ], // (498_H_CA, 498_H_C, 499_V_N)
     [  13.39499, 123.60153,  14.80645, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "H", 498, "CNCA" ], // (498_H_C, 499_V_N, 499_V_CA)
     [  14.80645, 110.85704,  15.30921, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "V", 499, "NCAC" ], // (499_V_N, 499_V_CA, 499_V_C)
     [  15.30921, 121.15316,  12.31572, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "V", 499, "CACO" ], // (499_V_CA, 499_V_C, 499_V_O)
     [  15.30921, 116.64412,  13.36851, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "V", 499, "CACN" ], // (499_V_CA, 499_V_C, 500_L_N)
     [  13.36851, 124.32061,  14.66764, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "V", 499, "CNCA" ], // (499_V_C, 500_L_N, 500_L_CA)
     [  14.66764, 113.54123,  15.38251, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "L", 500, "NCAC" ], // (500_L_N, 500_L_CA, 500_L_C)
     [  15.38251, 119.53917,  12.33013, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "L", 500, "CACO" ], // (500_L_CA, 500_L_C, 500_L_O)
     [  15.38251, 118.33943,  13.38837, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "L", 500, "CACN" ], // (500_L_CA, 500_L_C, 501_R_N)
     [  13.38837, 124.39461,  14.71909, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "L", 500, "CNCA" ], // (500_L_C, 501_R_N, 501_R_CA)
     [  14.71909, 111.49998,  15.33168, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "R", 501, "NCAC" ], // (501_R_N, 501_R_CA, 501_R_C)
     [  15.33168, 118.49783,  12.27333, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "R", 501, "CACO" ], // (501_R_CA, 501_R_C, 501_R_O)
     [  15.33168, 118.63538,  13.36515, "Csb", "Cdb", "Nsb", 0, 0, 1, 0, StdBond, "R", 501, "CACN" ], // (501_R_CA, 501_R_C, 502_G_N)
     [  13.36515, 123.89593,  14.67625, "Cdb", "Nsb", "Csb", 0, 0, 1, 0, StdBond, "R", 501, "CNCA" ], // (501_R_C, 502_G_N, 502_G_CA)
     [  14.67625, 112.46147,  15.30232, "Nsb", "Csb", "Cdb", 0, 0, 1, 0, StdBond, "G", 502, "NCAC" ], // (502_G_N, 502_G_CA, 502_G_C)
     [  15.30232, 119.84700,  12.29579, "Csb", "Cdb", "Odb", 0, 0, 1, 0, StdBond, "G", 502, "CACO" ], // (502_G_CA, 502_G_C, 502_G_O)
   ],

[  // chain - world transform for each residue
     [ 0, "260L", //(260_L_N, 260_L_CA, 260_L_C)
[ [ 1.0, 0.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0, 0.0 ], [ 0.0, 0.0, 1.0, 0.0 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 1, "261D", //(261_D_N, 261_D_CA, 261_D_C)
[ [ 0.020994799902191708, -0.9572037972214138, 0.28865222840292315, 9.89873760614736 ], [ 0.39981511112032425, 0.27266299885026557, 0.8751015746630918, -8.660978146947828 ], [ -0.9163553324431029, 0.09703494031989943, 0.3884290476524687, 36.222969827610726 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 2, "262A", //(262_A_N, 262_A_CA, 262_A_C)
[ [ -0.6205222549285349, -0.2552431588904417, -0.741487060560085, 23.61948693105484 ], [ -0.7220972486488164, -0.18277046604844643, 0.6672110012838839, 27.450183169302917 ], [ -0.30582297924161245, 0.9494450413689369, -0.07089724104449946, 37.57061028340331 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 3, "263L", //(263_L_N, 263_L_CA, 263_L_C)
[ [ 0.5395659288289276, 0.4595627628986497, -0.7054577771943096, -9.843041232683074 ], [ -0.7331078174008117, 0.6684820249533752, -0.1252386137740806, 44.181076142064406 ], [ 0.4140308400503456, 0.5847511002736074, 0.6975991787667261, 27.30549840216171 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 4, "264S", //(264_S_N, 264_S_CA, 264_S_C)
[ [ 0.5363821422154388, -0.005606885012855395, -0.8439566104681149, -35.01682547729275 ], [ -0.03492104472360057, -0.9992690050363078, -0.015555584500742542, 53.21856390011092 ], [ -0.8432524640626846, 0.037815584276856946, -0.5361858478531588, 54.28058454465295 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 5, "265P", //(265_P_N, 265_P_CA, 265_P_C)
[ [ 0.8807405388832337, 0.4063693520738333, -0.24322839649753894, -72.45720508680996 ], [ 0.354858488282156, -0.2261262206094736, 0.9071617196768053, 48.7420705205479 ], [ 0.31364240220635675, -0.8852857629308, -0.34336214335633647, 45.543313557305765 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 6, "266E", //(266_E_N, 266_E_CA, 266_E_C)
[ [ 0.36459663601006265, -0.7814565650309722, 0.506354549678523, -74.10561008474731 ], [ -0.9301179852317943, -0.2798471646833636, 0.23783628395816747, 87.21227825192499 ], [ -0.044156840449624445, -0.5576837826122222, -0.8288781406532749, 43.187819537016225 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 7, "267Q", //(267_Q_N, 267_Q_CA, 267_Q_C)
[ [ -0.7462146180836782, -0.5801763072882592, -0.326434061059182, -44.98307148652541 ], [ -0.3800618519891879, 0.7738786483800888, -0.5066209887519026, 89.9406915702321 ], [ 0.5465498444065088, -0.25398285383617214, -0.797982441872285, 18.23825738037941 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 8, "268L", //(268_L_N, 268_L_CA, 268_L_C)
[ [ 0.10103147845552131, 0.574489558867066, -0.8122526621155608, -60.154119687853395 ], [ 0.7495779146071481, 0.4928676064763088, 0.44183081877514313, 61.84976214759209 ], [ 0.654160217602913, -0.6534854774507528, -0.3808295425340202, -3.2144024971763603 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 9, "269V", //(269_V_N, 269_V_CA, 269_V_C)
[ [ 0.9675369207710052, 0.031220506126947864, 0.2507939132876609, -93.50314445162948 ], [ -0.2027629350352431, -0.4964677833691052, 0.8440420204305358, 81.45386250644017 ], [ 0.15086251728264158, -0.8674935273958789, -0.4740205489273648, -5.421588263040082 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 10, "270L", //(270_L_N, 270_L_CA, 270_L_C)
[ [ -0.22202253673629296, -0.8931094428496614, 0.39123077112357696, -74.06681111722123 ], [ -0.9622846262990473, 0.13601186934004914, -0.23560362770366738, 113.28653711233797 ], [ 0.15720779614791483, -0.42878467147438076, -0.8896231867138598, -15.135380719345408 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 11, "271T", //(271_T_N, 271_T_CA, 271_T_C)
[ [ -0.6778176389958055, -0.05081154624466498, -0.7334721774098716, -54.87625300809844 ], [ 0.2453856711829139, 0.9247722316497978, -0.2908301771611119, 94.0276870195002 ], [ 0.6930722333525292, -0.3771133865797914, -0.6143585052866593, -42.53025956797498 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 12, "272L", //(272_L_N, 272_L_CA, 272_L_C)
[ [ 0.6667801377297892, 0.5746331666076715, -0.47455344458077997, -88.49197423589084 ], [ 0.5870697905173763, -0.012748769019662931, 0.8094359331969191, 79.7663859733482 ], [ 0.45907876120653546, -0.81831179431982, -0.34585039870751844, -55.100569237797274 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 13, "273L", //(273_L_N, 273_L_CA, 273_L_C)
[ [ 0.6534951560573843, -0.5226024214360536, 0.5475680689363728, -102.23977972101085 ], [ -0.7564229958611498, -0.47737663476592834, 0.4471418118583862, 115.86013192194483 ], [ 0.027718808451626612, -0.7063980872629334, -0.7072718077014601, -56.89885023256505 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 14, "274E", //(274_E_N, 274_E_CA, 274_E_C)
[ [ -0.6737482577126348, -0.7371566099453926, 0.05160830982511863, -71.10011732351283 ], [ -0.6213281036904617, 0.527307547001544, -0.5795671992440818, 128.33606692565303 ], [ 0.40001834057153424, -0.42254808399640237, -0.8132886596512737, -76.23496298635028 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 15, "275A", //(275_A_N, 275_A_CA, 275_A_C)
[ [ -0.32768902517018456, 0.5859241897861622, -0.7411563577319185, -73.59260660312236 ], [ 0.7723038463809454, 0.6179889292045104, 0.14709334534865826, 97.14519842693002 ], [ 0.5442119730842121, -0.5241970308997295, -0.6550196952364844, -98.3142928435321 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 16, "276E", //(276_E_N, 276_E_CA, 276_E_C)
[ [ 0.9760855908876918, 0.17836945887421085, 0.12426284803732335, -108.72557572600451 ], [ 0.1319703355001361, -0.9404400534058932, 0.31329911665675414, 105.2040954298169 ], [ 0.17274475334842224, -0.28940774366081173, -0.9414894625537509, -111.96373244941972 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 17, "277P", //(277_P_N, 277_P_CA, 277_P_C)
[ [ -0.274174094500711, 0.8202990162685205, -0.5019343481109981, -112.43963131325695 ], [ 0.033221847670896316, 0.5297023870024475, 0.8475327073577993, 107.56919532016292 ], [ 0.9611060684139119, 0.2156963261458265, -0.17248252127432015, -150.0543615906769 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 18, "278P", //(278_P_N, 278_P_CA, 278_P_C)
[ [ 0.1635523068852371, -0.45021139543609384, -0.8778156653489401, -122.5806596171958 ], [ -0.9804249753233318, 0.024698961037434165, -0.19533773083026076, 141.07425173330031 ], [ 0.10962440729491461, 0.8925803385371854, -0.43735892420537414, -166.21743015020792 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 19, "279H", //(279_H_N, 279_H_CA, 279_H_C)
[ [ 0.980308782497815, 0.03822978552865986, 0.19373480445207653, -157.24517116737266 ], [ 0.04765089965813519, -0.9978855107536115, -0.044202931913778123, 147.15804456645648 ], [ 0.19163528568460975, 0.05256416009446398, -0.9800576137932614, -181.07010285703282 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 20, "280V", //(280_V_N, 280_V_CA, 280_V_C)
[ [ -0.36407745225748617, 0.5654422621742401, -0.7400828716468013, -159.60674170813718 ], [ 0.33473898114354117, 0.8209644493874129, 0.46256587351965983, 134.92311425070295 ], [ 0.8691360211501316, -0.07932478168454009, -0.4881702118626435, -217.46626280694227 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 21, "281L", //(281_L_N, 281_L_CA, 281_L_C)
[ [ 0.6890737738381947, -0.7236808570139023, -0.038253776284342914, -177.61987447977225 ], [ -0.6299178936095552, -0.5720247085269878, -0.5253486272420591, 153.76671123293104 ], [ 0.35830263956450253, 0.38610069933404506, -0.8500267457291397, -245.9515780350317 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 22, "282I", //(282_I_N, 282_I_CA, 282_I_C)
[ [ 0.07191281264224415, 0.9651336687792459, -0.2516854162771345, -192.218362969096 ], [ 0.8521903115684879, 0.07166406571331609, 0.5183010076724152, 136.5068069061243 ], [ 0.518266553278026, -0.2517563565714886, -0.8173239973714002, -277.46916897939553 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 23, "283S", //(283_S_N, 283_S_CA, 283_S_C)
[ [ 0.3384375524605949, -0.3423740872770737, -0.8764930161990298, -211.47464234724015 ], [ -0.7294180174555152, 0.49300972997399767, -0.4742264880436684, 163.54644213610274 ], [ 0.5944824462469415, 0.7998258501149837, -0.08288082161812604, -296.40302801921143 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 24, "284R", //(284_R_N, 284_R_CA, 284_R_C)
[ [ 0.8790042089654924, 0.3609456564390544, -0.3115603211558383, -249.3039264032221 ], [ 0.3729755761945046, -0.9275724937805497, -0.0223268525347699, 156.44357401767354 ], [ -0.29705356450197035, -0.09657899295142909, -0.9499640403389812, -292.3572781115818 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 25, "285P", //(285_P_N, 285_P_CA, 285_P_C)
[ [ 0.1992099574001478, 0.6264150488678728, -0.7536043918558939, -267.08126822556943 ], [ 0.27409756038076083, 0.7027060608878989, 0.6565628068850136, 145.03308387906057 ], [ 0.9408431964285794, -0.33735497409001247, -0.03171279222617503, -324.642877411658 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 26, "286S", //(286_S_N, 286_S_CA, 286_S_C)
[ [ 0.6184936251173931, 0.002798064973949348, -0.7857848347490216, -294.73211568760547 ], [ -0.7436590895014654, 0.3251253953318016, -0.5841785993274584, 166.59898206226168 ], [ 0.25384403536620115, 0.9456667743675577, 0.20316879082738465, -339.6380983789283 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 27, "287A", //(287_A_N, 287_A_CA, 287_A_C)
[ [ 0.7280003384479335, 0.390606671211459, -0.5634198573220482, -329.1232157127756 ], [ 0.5030334424508978, -0.8627085746628368, 0.05187746118707099, 156.1182844683301 ], [ -0.4658034596218596, -0.3211858396759108, -0.8245403528002657, -326.648288078582 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 28, "288P", //(288_P_N, 288_P_CA, 288_P_C)
[ [ 0.4880747310618513, 0.5157119684319177, 0.7041478697794791, -349.6489945651448 ], [ 0.5912464238477042, 0.39810322806682374, -0.7013854048168349, 135.3589306865982 ], [ -0.6420363877430953, 0.7586534026937802, -0.11060873109726164, -333.7814606301792 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 29, "289F", //(289_F_N, 289_F_CA, 289_F_C)
[ [ -0.8415799255792766, 0.38191247599925227, -0.3819503757533166, -326.84216216662844 ], [ 0.5170151343488937, 0.36492385618820683, -0.7742905979274993, 105.77122871292391 ], [ -0.15632843540499788, -0.8491015486152699, -0.5045671218212457, -325.44550586900436 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 30, "290T", //(290_T_N, 290_T_CA, 290_T_C)
[ [ 0.3319327955880142, -0.7792353505454368, 0.5316134758200319, -326.44447760747397 ], [ 0.7365552031749271, -0.13798481022717338, -0.6621530222105105, 74.56000753806119 ], [ 0.5893276269521627, 0.6113529754625401, 0.5281481681346925, -347.41576266232676 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 31, "291E", //(291_E_N, 291_E_CA, 291_E_C)
[ [ -0.1217685215133152, 0.11950031044112641, 0.9853385727621423, -314.92512015203766 ], [ 0.8263974682300202, 0.5620621930823937, 0.03396050073869158, 40.527809618789675 ], [ -0.5497632687543558, 0.8184166218448183, -0.1671962362514229, -333.2380523654309 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 32, "292A", //(292_A_N, 292_A_CA, 292_A_C)
[ [ -0.8459926368016778, 0.28368587055904937, 0.4514629390354194, -281.9377423243744 ], [ 0.2709303697063457, -0.5005464023260604, 0.822222618205812, 43.082748365830895 ], [ 0.45923108915678534, 0.8179093017954937, 0.34659945295438754, -352.77523542232706 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 33, "293S", //(293_S_N, 293_S_CA, 293_S_C)
[ [ -0.3801602610661274, 0.9133046018759546, 0.14612624712329364, -275.5071388411766 ], [ -0.9240998214846868, -0.3817071583047383, -0.018416439147982782, 79.92852933788599 ], [ 0.0389576159191472, -0.14203643719526848, 0.9890945125066516, -342.54819306206105 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 34, "294M", //(294_M_N, 294_M_CA, 294_M_C)
[ [ 0.21054696777184456, 0.6367318361336888, 0.7417833532884808, -281.8396227294378 ], [ -0.0898271486926042, 0.7681835397593103, -0.633896783869908, 71.61258962509874 ], [ -0.9734480251764028, 0.0668327621505905, 0.21893452031026073, -305.6227372085101 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 35, "295M", //(295_M_N, 295_M_CA, 295_M_C)
[ [ -0.34958798380621153, 0.162798567699243, 0.9226510000717405, -261.2569719120013 ], [ 0.8786031187662339, 0.39893089069559223, 0.26250848394493576, 39.396591850645706 ], [ -0.32533798006466313, 0.9024138578301513, -0.282496775067519, -304.9777727713571 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 36, "296M", //(296_M_N, 296_M_CA, 296_M_C)
[ [ -0.8081008584546073, 0.4661153278323246, 0.36015205639381176, -232.3760803871734 ], [ -0.10821633490999166, -0.7184890247942536, 0.6870682252213757, 55.97959821497849 ], [ 0.5790183308182698, 0.5162460870650423, 0.6310528901498332, -324.77461243799445 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 37, "297S", //(297_S_N, 297_S_CA, 297_S_C)
[ [ -0.20088721376687796, 0.9431621299908775, 0.2647442613090072, -231.36346445652302 ], [ -0.8851703355576892, -0.058989472094720566, -0.46151242586811686, 84.02354899907836 ], [ -0.4196639183843049, -0.32705571197130523, 0.8467093697800173, -298.4851061658997 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 38, "298L", //(298_L_N, 298_L_CA, 298_L_C)
[ [ 0.14070947443703768, 0.545396624170611, 0.82628286085756, -232.03727371946246 ], [ 0.3619751418254863, 0.7484654174283708, -0.5556739292194602, 58.5713957248714 ], [ -0.9215068315015748, 0.37728244228571783, -0.0921038448646543, -269.7468616325062 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 39, "299T", //(299_T_N, 299_T_CA, 299_T_C)
[ [ -0.5922403788252151, 0.12894529904878402, 0.795376919166118, -204.3511729885629 ], [ 0.7923039030457846, -0.08646185372534688, 0.6039692650034377, 36.026584708225 ], [ 0.14664876033360053, 0.9878752237522695, -0.0509576627115277, -284.92512921912305 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 40, "300K", //(300_K_N, 300_K_CA, 300_K_C)
[ [ -0.6922520373101133, 0.6954515912851561, 0.19270236381254258, -183.00087525556756 ], [ -0.5870904033859897, -0.698002137643898, 0.410021797099404, 67.28335314254558 ], [ 0.4196569731245668, 0.17070471587947703, 0.8914863570938333, -292.67877816085416 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 41, "301L", //(301_L_N, 301_L_CA, 301_L_C)
[ [ 0.016691171348163725, 0.8693847961945911, 0.493853703990074, -187.68860261703125 ], [ -0.5384902835184854, 0.4239873280909381, -0.7281915683214745, 78.370601818985 ], [ -0.8424663906383492, -0.25378105083712726, 0.4752321105321019, -256.20556052941214 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 42, "302A", //(302_A_N, 302_A_CA, 302_A_C)
[ [ -0.025766507440532778, 0.27987742108352387, 0.9596899063040899, -179.59351123827733 ], [ 0.7767762297884581, 0.609890881704488, -0.15700892092282825, 43.52642097182628 ], [ -0.629249374993687, 0.7414187356557415, -0.2331168858935281, -242.45261391307784 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 43, "303D", //(303_D_N, 303_D_CA, 303_D_C)
[ [ -0.7807420690857604, 0.20325861513040655, 0.5908703385134129, -147.2332971051634 ], [ 0.36249925601085087, -0.6228866919021308, 0.6932578585510588, 41.13321988848485 ], [ 0.5089559027570714, 0.7554456330051087, 0.41263274787908577, -263.3956427465701 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 44, "304K", //(304_K_N, 304_K_CA, 304_K_C)
[ [ -0.5115933800020751, 0.8572629522408997, 0.05807274966169093, -135.21473559495027 ], [ -0.8542503270694078, -0.5002024492115851, -0.1416117527062683, 75.32581411651307 ], [ -0.09235037758375474, -0.12205630060733715, 0.9882174190137458, -249.94765431287078 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 45, "305E", //(305_E_N, 305_E_CA, 305_E_C)
[ [ 0.22753912629113152, 0.7988974987504124, 0.5567663176746627, -145.40273596108042 ], [ -0.12512455170694167, 0.5910115627841511, -0.7968997297123223, 66.74671223030263 ], [ -0.9656965323365334, 0.11166073233581382, 0.23443994601480433, -213.88810769136552 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 46, "306L", //(306_L_N, 306_L_CA, 306_L_C)
[ [ -0.27250172445075593, -0.06280286409684903, 0.9601034373610986, -126.6961431459039 ], [ 0.9557518544011699, 0.09726524112461683, 0.27762900726994455, 32.94192860891279 ], [ -0.11082058915241814, 0.9932750239134189, 0.033519007886691446, -217.31979654907335 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 47, "307V", //(307_V_N, 307_V_CA, 307_V_C)
[ [ -0.9260974470090513, 0.3312158578845113, 0.18066425802893202, -93.69280430768369 ], [ -0.21101501874806514, -0.8516731320541272, 0.4797140168890908, 51.38420153350517 ], [ 0.3127557841288807, 0.40613905454036353, 0.8586238337426865, -225.64569467731792 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 48, "308H", //(308_H_N, 308_H_CA, 308_H_C)
[ [ -0.07473851802634662, 0.9963978798613604, -0.04006767938141127, -99.91239309310936 ], [ -0.7731441540603288, -0.08327574731356183, -0.6287394269104702, 73.24571153623572 ], [ -0.6298112979024465, -0.016012960911521337, 0.7765831018746686, -194.7597716402301 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 49, "309M", //(309_M_N, 309_M_CA, 309_M_C)
[ [ 0.3994154884853496, 0.30151384284991706, 0.8657694093278447, -109.04112679155585 ], [ 0.5764421895771826, 0.6517318837480648, -0.49290968116041556, 42.199502703307814 ], [ -0.7128686201772899, 0.6959417750617322, 0.08650650894847361, -173.304569468401 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 50, "310I", //(310_I_N, 310_I_CA, 310_I_C)
[ [ -0.7022987561490909, -0.20243961112839787, 0.682491509806256, -76.23802110718634 ], [ 0.696508946752586, -0.393566539353822, 0.5999838882792268, 24.249088260933792 ], [ 0.14714531662626246, 0.8967293811107337, 0.41740229137813495, -183.02807422693093 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 51, "311S", //(311_S_N, 311_S_CA, 311_S_C)
[ [ -0.7171421388956705, 0.6748792873475506, -0.17391118460123617, -57.92511089613267 ], [ -0.6962570154261501, -0.7047262010486731, 0.13633469853773464, 56.691992884316065 ], [ -0.03055030425401415, 0.21885823965473036, 0.9752782935373984, -173.56368842163002 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 52, "312W", //(312_W_N, 312_W_CA, 312_W_C)
[ [ 0.4453091335576187, 0.8401524163008838, 0.3095863255280333, -77.61936223164186 ], [ -0.36716734627306385, 0.48669782282139007, -0.7926622036477549, 59.241146489202805 ], [ -0.8166320563148124, 0.23930972952369753, 0.5252075189430705, -140.30868767079744 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 53, "313A", //(313_A_N, 313_A_CA, 313_A_C)
[ [ 0.08149178161940442, -0.1811497465117875, 0.9800733946328767, -70.5028973942027 ], [ 0.9112221823409573, 0.4119150380794904, 0.00036852917790088, 22.551234866982487 ], [ -0.4037737286380558, 0.8930345854123943, 0.19863535767590473, -130.84192728569468 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 54, "314K", //(314_K_N, 314_K_CA, 314_K_C)
[ [ -0.9633136831056955, -0.010904921200525833, 0.26815635482857736, -32.55111024329574 ], [ 0.19075081835621757, -0.7306797636364069, 0.6555312412911211, 28.954380885920592 ], [ 0.18878790543298282, 0.6826332585699696, 0.7059539369224965, -135.29223012505477 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 55, "315K", //(315_K_N, 315_K_CA, 315_K_C)
[ [ -0.28855719711581734, 0.9211459313237967, -0.2611989992290908, -34.9580991626042 ], [ -0.8436744177708316, -0.37360793499930367, -0.38552637744336565, 57.19550111527122 ], [ -0.452712072725824, 0.10912050270765479, 0.8849545158352021, -108.87296370560462 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 56, "316I", //(316_I_N, 316_I_CA, 316_I_C)
[ [ 0.5507815061239422, 0.239243415500975, 0.7996263631542359, -50.627691951140555 ], [ 0.3056393719259023, 0.8336813827535299, -0.45995643965369426, 33.98044819351959 ], [ -0.7766751617250174, 0.4977327999937219, 0.3860540803693541, -82.3826672902127 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 57, "317P", //(317_P_N, 317_P_CA, 317_P_C)
[ [ -0.8673209820486044, 0.49451946088680804, -0.05661110228979193, -23.06295434360844 ], [ 0.07155798875039536, 0.011326975341986128, -0.9973721240718533, 23.459132150136845 ], [ -0.49257869253982434, -0.8690927467387026, -0.04521094139624152, -57.4552725601077 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 58, "318G", //(318_G_N, 318_G_CA, 318_G_C)
[ [ -0.22441483385929603, 0.9166242512986575, -0.33081409322318794, -11.87389843250024 ], [ 0.9025711191210832, 0.06751224303833038, -0.42521461871441907, -12.987338892638293 ], [ -0.3674280300581597, -0.39400771435481935, -0.8424693251130697, -61.6939879832722 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 59, "319F", //(319_F_N, 319_F_CA, 319_F_C)
[ [ 0.6788574824722526, 0.09743731179666605, 0.727776400250293, -29.588474599287792 ], [ 0.33952003989403146, -0.9204886122026604, -0.19346022153294362, -18.031376356075533 ], [ 0.6510596447344826, 0.37842659139529244, -0.6579625019126468, -95.6344191764093 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 60, "320V", //(320_V_N, 320_V_CA, 320_V_C)
[ [ -0.5360079687719181, -0.5059654270606232, 0.6757917164573444, 2.6221998246740803 ], [ -0.16849570926576604, -0.7202752850602068, -0.6729135974926235, -13.634963887110766 ], [ 0.8272270869428956, -0.47455505513143487, 0.30082028900519975, -116.55166204500014 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 61, "321E", //(321_E_N, 321_E_CA, 321_E_C)
[ [ -0.7542399851062018, 0.5384013019695526, -0.3758271981968404, 20.45942570694524 ], [ 0.39103590587936005, -0.09148241131132014, -0.9158176066956645, -38.46368857598663 ], [ -0.5274589701390098, -0.8377081869351386, -0.14153454829736134, -92.93209235065812 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 62, "322L", //(322_L_N, 322_L_CA, 322_L_C)
[ [ 0.572526882139162, 0.7034428020349514, 0.421166467674275, -3.441213913567484 ], [ 0.8189802875659501, -0.5148099694204151, -0.25346002438993187, -66.26333585890927 ], [ 0.038526066583596995, 0.49003971231989685, -0.8708483349832099, -105.08247828801409 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 63, "323S", //(323_S_N, 323_S_CA, 323_S_C)
[ [ -0.5627135440014546, 0.6346033435015086, -0.5297471697083705, 7.0527429970515705 ], [ 0.456021490043453, 0.7728081717254939, 0.4413750449819842, -86.8209572437566 ], [ 0.6894910209827767, 0.006791622118900699, -0.7242623874350526, -135.76972579513665 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 64, "324L", //(324_L_N, 324_L_CA, 324_L_C)
[ [ 0.21141831689513113, -0.5714994879593736, -0.7929001390738961, -2.371430160781529 ], [ -0.6481290255797922, 0.5252485974820098, -0.5514006502028908, -71.93967211526598 ], [ 0.7315948752432722, 0.6304777919208221, -0.25935823181928425, -170.0681976202632 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 65, "325F", //(325_F_N, 325_F_CA, 325_F_C)
[ [ 0.5215367483095003, -0.40703812563543434, -0.7498795799606397, -24.636415828166854 ], [ 0.7500074537954758, 0.6377318865628209, 0.1754618480247274, -102.92344149040588 ], [ 0.4068024574827199, -0.6539250760906923, 0.6378822426167718, -174.81267794131824 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 66, "326D", //(326_D_N, 326_D_CA, 326_D_C)
[ [ 0.5821087394021268, -0.7067497026820659, -0.40206252408111637, -43.14552699598765 ], [ 0.11381035615961947, -0.4187873582184067, 0.9009241651922033, -95.10503580469447 ], [ -0.8051065881872334, -0.5701947091610544, -0.16334434579211254, -142.05007624457025 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 67, "327Q", //(327_Q_N, 327_Q_CA, 327_Q_C)
[ [ 0.08552737122489039, -0.7483235785750053, -0.6577969979561001, -45.915083953584066 ], [ -0.9931253333593145, -0.11698878697598253, 0.0039618129484455915, -57.48983466512137 ], [ -0.07991959091055262, 0.652936019431143, -0.7531848468457819, -150.47597107425298 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 68, "328V", //(328_V_N, 328_V_CA, 328_V_C)
[ [ 0.3557230468220761, -0.4622973066406911, -0.8122452303537288, -61.80200024434587 ], [ -0.08746493587809402, 0.8488096160545179, -0.5214135793065068, -65.42568176380503 ], [ 0.9304896554779424, 0.25652180407549136, 0.2615063385128461, -184.77429012101427 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 69, "329R", //(329_R_N, 329_R_CA, 329_R_C)
[ [ 0.54025823578017, -0.5247112419804834, -0.6578747230369905, -83.62428275476526 ], [ 0.8398851458700413, 0.28782892755190415, 0.4601602440576328, -93.53460148992463 ], [ -0.05209587717427992, -0.8011445693528848, 0.596199126616119, -169.91291768020065 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 70, "330L", //(330_L_N, 330_L_CA, 330_L_C)
[ [ 0.40251091927252647, -0.7543026575646223, -0.5186641115956778, -97.18483757301316 ], [ -0.36945903319260154, -0.6522538274921567, 0.6618647651251878, -69.32867997555155 ], [ -0.8375470032534137, -0.07478265382279542, -0.5412232183013272, -143.02087556029198 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 71, "331L", //(331_L_N, 331_L_CA, 331_L_C)
[ [ 0.1883334520253373, -0.6892992954934263, -0.6995691474618404, -104.68997209519783 ], [ -0.8890471238291313, 0.18299515476586314, -0.41965221903782074, -41.39355037507373 ], [ 0.41728374334422524, 0.7009844895320324, -0.5783554469149113, -168.55439489807156 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 72, "332E", //(332_E_N, 332_E_CA, 332_E_C)
[ [ 0.3297447885800482, -0.5122068801834994, -0.7930400281807883, -121.04694566661863 ], [ 0.4571937921329771, 0.8215861515686936, -0.34054372991674875, -66.76808312676292 ], [ 0.8259795462596924, -0.2502804575735583, 0.5050915577570099, -192.71169544339048 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 73, "333S", //(333_S_N, 333_S_CA, 333_S_C)
[ [ 0.48649607455264154, -0.5059864955549106, -0.7122494196283569, -143.27727638747845 ], [ 0.7444749761550639, -0.18657477292948416, 0.6410513739047976, -79.74397737841012 ], [ -0.4572511118890668, -0.8421208466855408, 0.28592638957228184, -163.843403066034 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 74, "334C", //(334_C_N, 334_C_CA, 334_C_C)
[ [ 0.5025795977208266, -0.8513645096541048, -0.15030708451765645, -158.29030298547684 ], [ -0.3464557256067792, -0.35762702795268436, 0.8672204673968495, -47.36914183786113 ], [ -0.7920746039036047, -0.38377256360917167, -0.47469615679062754, -150.4539558904468 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 75, "335W", //(335_W_N, 335_W_CA, 335_W_C)
[ [ -0.3763274559297251, 0.13650008644903427, -0.916376217670923, -156.10618115118433 ], [ -0.7948907502722838, 0.46053709531130355, 0.3950370602535832, -19.505702059333224 ], [ 0.47594783437371213, 0.8770822710793061, -0.06481009731006374, -176.87021292236702 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 76, "336M", //(336_M_N, 336_M_CA, 336_M_C)
[ [ 0.8803197944910096, 0.3034785122587582, -0.3646064344255558, -193.92841210943115 ], [ -0.2263454637141267, 0.944173809080766, 0.23938159766772843, -17.042229763678336 ], [ 0.41689901712926036, -0.1282053463905775, 0.8998658781582557, -182.70260070588023 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 77, "337E", //(337_E_N, 337_E_CA, 337_E_C)
[ [ 0.4135759625504516, -0.8262843041916894, 0.3823730793961877, -199.04474499325934 ], [ 0.09338018438806375, 0.4562541332000037, 0.884936329405435, -16.794744137781812 ], [ -0.9056682970956413, -0.3302823255706685, 0.26585432298819633, -144.39725261935078 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 78, "338V", //(338_V_N, 338_V_CA, 338_V_C)
[ [ -0.6424121957297734, -0.6544765173575237, -0.3987067330823728, -173.82255599459683 ], [ -0.6344759642714939, 0.1624169717887661, 0.7556857667289528, 12.193594355609843 ], [ -0.42982184860638417, 0.7384315916200594, -0.5195882629140658, -140.4866340131619 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 79, "339L", //(339_L_N, 339_L_CA, 339_L_C)
[ [ 0.14136812387039469, 0.36991491229451745, -0.9182472495006474, -190.08801248285513 ], [ -0.6600072024583155, 0.7265527588879788, 0.19108003887223496, 31.76378428852178 ], [ 0.7378384282866494, 0.5790371717035152, 0.3468579068279704, -169.47774368451985 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 80, "340M", //(340_M_N, 340_M_CA, 340_M_C)
[ [ 0.9975413376050412, -0.06586871215013747, 0.02392890571730367, -225.55493064572585 ], [ 0.04979740495826537, 0.9064777850545703, 0.4193068609765203, 26.035302811278193 ], [ -0.0493102243816412, -0.41708432955727637, 0.9075291531456108, -155.44751960253578 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 81, "341M", //(341_M_N, 341_M_CA, 341_M_C)
[ [ -0.08794018654410492, -0.9775402161607012, 0.1915245398873932, -213.2443306815099 ], [ -0.08833410383686678, 0.19916479939885418, 0.9759766743010572, 35.2149794854654 ], [ -0.9922013957306882, 0.06890942220700898, -0.10386472857120704, -120.10677636594096 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 82, "342G", //(342_G_N, 342_G_CA, 342_G_C)
[ [ -0.5435161144918148, -0.20873903127045068, -0.813030288557563, -199.33653509288322 ], [ -0.828852295678745, 0.2865188652243842, 0.48053180104817783, 67.88318476879334 ], [ 0.132642773025125, 0.935058798522748, -0.3287414456216316, -135.17210300863889 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 83, "343L", //(343_L_N, 343_L_CA, 343_L_C)
[ [ 0.6579121529015182, 0.36743074438812745, -0.6573783135629558, -232.1620702042151 ], [ -0.3129123988588504, 0.9273581561212011, 0.20516500899492035, 74.39693710080857 ], [ 0.6850090727172079, 0.07072127228686267, 0.7250938366454637, -154.33244236633902 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 84, "344M", //(344_M_N, 344_M_CA, 344_M_C)
[ [ 0.8191507035398627, -0.490054731815132, 0.2980578546452802, -254.0464202078293 ], [ 0.15495276355356832, 0.6893945959113441, 0.7076190586716461, 69.89422842863846 ], [ -0.5522515422859947, -0.5334617614729604, 0.6406534032446144, -123.01629005363806 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 85, "345W", //(345_W_N, 345_W_CA, 345_W_C)
[ [ -0.470886544094323, -0.8717181318767839, -0.13554837198636974, -230.78914770270677 ], [ -0.42881180276457176, 0.09189125598974615, 0.8987082034133353, 92.76748618585852 ], [ -0.7709645260326206, 0.4813143418078013, -0.4170733795981467, -102.68598956891006 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 86, "346R", //(346_R_N, 346_R_CA, 346_R_C)
[ [ -0.20396193849128688, 0.11102422165291014, -0.9726629168695813, -232.89092641536243 ], [ -0.8097637345929158, 0.5392206729639262, 0.23135202611272782, 120.39757127984183 ], [ 0.5501656312284604, 0.8348141637841898, -0.02007710537636645, -129.686050902843 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 87, "347S", //(347_S_N, 347_S_CA, 347_S_C)
[ [ 0.879255152959383, 0.10219818869658989, -0.4652589668362217, -271.1152031742949 ], [ -0.06515510167288953, 0.9933366509342318, 0.09506371882458318, 115.43518605370619 ], [ 0.47187412380884763, -0.0532712693475591, 0.8800551023326506, -129.53966856131802 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 88, "348I", //(348_I_N, 348_I_CA, 348_I_C)
[ [ 0.1727445917198085, -0.7370295493321268, -0.6534116232841582, -276.8744801553077 ], [ 0.34890043839407575, -0.5745830583917663, 0.7403531543105533, 116.48369173057351 ], [ -0.9211014005634659, -0.3558676051854764, 0.15789381704036257, -91.4295284411156 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 89, "349D", //(349_D_N, 349_D_CA, 349_D_C)
[ [ 0.3006381831308608, -0.6775446580026072, -0.6712301537147338, -291.5045933882474 ], [ -0.8356811123558509, -0.5263072784131524, 0.1569640950696311, 152.3076537119095 ], [ -0.4596234995030854, 0.5137449611446039, -0.7244393374210616, -91.06446013270443 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 90, "350H", //(350_H_N, 350_H_CA, 350_H_C)
[ [ 0.3770914953135705, -0.6153415956390766, -0.6922114740731553, -312.5612796789909 ], [ -0.38732187686476793, 0.5741214403307443, -0.7213642183075795, 146.40946194266436 ], [ 0.8412988576386687, 0.5401289590726522, -0.021838949309987046, -123.0440467245993 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 91, "351P", //(351_P_N, 351_P_CA, 351_P_C)
[ [ 0.6226021275312016, 0.6607547781948133, -0.4192489879371739, -344.5256503160493 ], [ 0.6822686973979403, -0.7207306739714949, -0.12270582768358017, 125.97947066924219 ], [ -0.3832440675921015, -0.20964355150902642, -0.8995407528100897, -114.88558590052166 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 92, "352G", //(352_G_N, 352_G_CA, 352_G_C)
[ [ 0.7442869483418749, 0.33497366772569476, 0.5777798719740374, -358.5720840880204 ], [ 0.33466790960614073, -0.9357300807508497, 0.11138494628029477, 107.79311155636145 ], [ 0.5779570302434215, 0.11046202020787044, -0.8085566234246064, -145.78147830846044 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 93, "353K", //(353_K_N, 353_K_CA, 353_K_C)
[ [ -0.34591812415534784, -0.3294748402718502, 0.8785140755893921, -327.7303562247541 ], [ -0.32389097837997, -0.8368291291872301, -0.441374719108161, 116.22336968206966 ], [ 0.8805880339323894, -0.4372222983467117, 0.18276043423976598, -167.79486729772884 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 94, "354L", //(354_L_N, 354_L_CA, 354_L_C)
[ [ -0.7557222075162628, 0.16120452820453812, 0.6347417153087537, -293.46403428180366 ], [ 0.6216965498412437, 0.481260966035864, 0.6179654379378519, 100.89041516991824 ], [ -0.20585758422326997, 0.8616269393749125, -0.46391989864757915, -175.1230398025657 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 95, "355I", //(355_I_N, 355_I_CA, 355_I_C)
[ [ -0.49470586609700035, -0.4184736598257186, 0.7616730939722668, -262.70522701479746 ], [ -0.8451492631770566, 0.02748425618717118, -0.533823321533554, 123.54813825750844 ], [ 0.2024569806168586, -0.9078129827744846, -0.36726932802181866, -179.6649256894465 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 96, "356F", //(356_F_N, 356_F_CA, 356_F_C)
[ [ -0.4381136288653914, 0.3523512251360543, 0.826985527290242, -242.3929255970119 ], [ 0.7630408472248623, -0.3405724554036836, 0.5493442163950306, 103.09199032321031 ], [ 0.4752105992807425, 0.8716989255272874, -0.11964894302002863, -204.79426601520157 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 97, "357A", //(357_A_N, 357_A_CA, 357_A_C)
[ [ -0.9043050302518625, 0.1808033551145654, -0.386707329954443, -220.72091885321188 ], [ -0.41691125230999115, -0.17936723079720565, 0.8910737366870696, 134.04125183535834 ], [ 0.09174649834458526, 0.9670250996202002, 0.23758164227492137, -211.05022906669672 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 98, "358P", //(358_P_N, 358_P_CA, 358_P_C)
[ [ 0.1456658478126668, -0.8385257097622381, -0.5250296133064799, -223.65821415864198 ], [ -0.848634040489761, -0.3786839249211548, 0.36934909006024647, 171.60171873543155 ], [ -0.5085289825595389, 0.391756453774592, -0.7667627761067409, -201.62469268274583 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 99, "359D", //(359_D_N, 359_D_CA, 359_D_C)
[ [ 0.23075890391006673, -0.6970637952265764, -0.6788611004104934, -235.87609452167263 ], [ -0.5539378770092573, 0.4794768284821815, -0.6806282387347335, 176.19621918799834 ], [ 0.7999394706355608, 0.5331079030862234, -0.2754864914770022, -238.02157004564987 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 100, "360L", //(360_L_N, 360_L_CA, 360_L_C)
[ [ 0.4921153152656172, -0.6955794770031553, -0.5234421722148868, -255.7193613586578 ], [ 0.6955202483203412, 0.6757692819228813, -0.24410543170082108, 142.76313246964773 ], [ 0.5235208693618938, -0.24393660812073778, 0.8163460237920661, -239.47332427550919 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 101, "361V", //(361_V_N, 361_V_CA, 361_V_C)
[ [ 0.5153553255051628, 0.4501555471689867, -0.7292248431221044, -285.3592033305761 ], [ 0.08483721472157994, -0.8735464308405421, -0.4792903923135061, 142.75224797923366 ], [ -0.8527669877942933, 0.1851394515459637, -0.4883767480229864, -215.43472937946973 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 102, "362L", //(362_L_N, 362_L_CA, 362_L_C)
[ [ 0.5436122368877763, -0.4456472735097831, -0.7112553996414823, -311.43853250192785 ], [ 0.5400297091568431, 0.8344172559419983, -0.1100715958556129, 115.86223389929862 ], [ 0.6425368854266996, -0.3242627801637469, 0.694262198500418, -223.0165719127264 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 103, "363D", //(363_D_N, 363_D_CA, 363_D_C)
[ [ 0.6986100473311339, 0.6186247825936569, -0.35950991659332804, -346.33697310027236 ], [ 0.0035174426586993025, -0.5054211753506569, -0.8628656112653368, 110.77727839354898 ], [ -0.7154939757713784, 0.6015420300096339, -0.35526828843390323, -208.68367974604593 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 104, "364R", //(364_R_N, 364_R_CA, 364_R_C)
[ [ 0.2460997058583777, 0.8426909591955061, -0.4788599816925382, -368.02683066938823 ], [ 0.9611467952710547, -0.2759093216822567, 0.008419272475592316, 79.12101878016938 ], [ -0.12502708793138084, -0.4623267172671238, -0.8778509177442396, -210.55499808300618 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 105, "365D", //(365_D_N, 365_D_CA, 365_D_C)
[ [ 0.7536794675019942, -0.04094675682141938, 0.6559654132435018, -386.91482486039183 ], [ -0.08050115669000432, -0.9962938081389668, 0.030302007120436906, 92.08054653916642 ], [ 0.652293510651058, -0.0756439751055111, -0.754182448080539, -241.50572953078765 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 106, "366E", //(366_E_N, 366_E_CA, 366_E_C)
[ [ -0.6647679660723055, -0.634216856621377, 0.394781623256947, -354.46515033513367 ], [ -0.3593450720559164, -0.19183212247024434, -0.9132751808615459, 99.52745029593017 ], [ 0.6549463111379576, -0.748978915301266, -0.10037885215980777, -260.89103837412387 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 107, "367G", //(367_G_N, 367_G_CA, 367_G_C)
[ [ -0.6051394405980279, 0.5559641185346914, -0.5698334461354908, -343.47990005287295 ], [ 0.7594298142360653, 0.1883334020673785, -0.622733399550002, 63.677710424625744 ], [ -0.23889875404045346, -0.8095890491895692, -0.5361836968336955, -251.75462662089038 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 108, "368K", //(368_K_N, 368_K_CA, 368_K_C)
[ [ 0.8148394239263457, 0.5055494811939801, 0.28364843606104834, -374.1303634485586 ], [ 0.5658001779299955, -0.8000646425817697, -0.1994159631145174, 48.65686385619394 ], [ 0.12612244792171842, 0.32298032409885685, -0.9379641988771378, -270.4887359277536 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 109, "369C", //(369_C_N, 369_C_CA, 369_C_C)
[ [ -0.09013035809896039, -0.7977192613489368, 0.5962553971427512, -355.5409083703948 ], [ -0.2134955777984225, -0.569298489256636, -0.7939262361143067, 51.23182706757037 ], [ 0.9727775474431299, -0.19885474649356175, -0.11899845792072429, -304.06998502144654 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 110, "370V", //(370_V_N, 370_V_CA, 370_V_C)
[ [ -0.7936024066689364, -0.3513165086393103, -0.4967614426329376, -334.0771379863679 ], [ 0.582182352750661, -0.201210933335613, -0.7877676487720959, 19.79187975773225 ], [ 0.17680194646829225, -0.9143800473893209, 0.364211752503027, -295.77742534745556 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 111, "371E", //(371_E_N, 371_E_CA, 371_E_C)
[ [ 0.7041661936424903, -0.5865655956526028, 0.4001134510707609, -349.55220733837416 ], [ 0.709670463117359, 0.5633588158329492, -0.42307762692216483, -15.25833901838236 ], [ 0.022755340248830743, 0.58186566028595, 0.8129665109155206, -293.5751728708143 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 112, "372G", //(372_G_N, 372_G_CA, 372_G_C)
[ [ -0.03828857062606414, -0.8394180850766502, 0.5421358351978427, -347.4935033560787 ], [ 0.7889641451844287, 0.30755562885128773, 0.5319258527043901, -30.66777038303974 ], [ -0.613245108396974, 0.44809241636981373, 0.6504949065281334, -258.305438136155 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 113, "373I", //(373_I_N, 373_I_CA, 373_I_C)
[ [ -0.8100243604972726, 0.29612615289842736, -0.5061322326927636, -326.7754545578186 ], [ -0.4790377259821984, -0.8319754365048119, 0.2798923545551608, -2.085065672802884 ], [ -0.33820613904361463, 0.46917605930205436, 0.815775969792569, -242.25573139830627 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 114, "374L", //(374_L_N, 374_L_CA, 374_L_C)
[ [ 0.6658955161562182, 0.7432176736459597, 0.06488953030598518, -356.5286097899972 ], [ -0.6117308522735638, 0.5937292149714491, -0.5227532722671814, 7.9897406439008325 ], [ -0.4270462807936662, 0.3084041323809953, 0.8500108029846578, -219.76237702380305 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 115, "375E", //(375_E_N, 375_E_CA, 375_E_C)
[ [ 0.1346199288876637, -0.5605800984595055, 0.8170847128403587, -351.1280008143515 ], [ 0.7088743287145273, 0.630648345575437, 0.3158794870080872, -21.43535251813941 ], [ -0.6923688762760849, 0.5366867032395112, 0.4822786764206938, -195.29137856567596 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 116, "376I", //(376_I_N, 376_I_CA, 376_I_C)
[ [ -0.9636856314509664, -0.24701863073182556, -0.10144850815229352, -314.80034108846013 ], [ -0.007761324848002271, -0.35383151600807144, 0.9352769750806622, -9.717886278477678 ], [ -0.26692651717571214, 0.902100357139273, 0.33906515609487264, -189.78438426864903 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 117, "377F", //(377_F_N, 377_F_CA, 377_F_C)
[ [ 0.004301677873234298, 0.8693060576945091, -0.4942554740446541, -326.338834706755 ], [ -0.9986656503235972, 0.02917076956125092, 0.042614376294286535, 26.977719319595373 ], [ 0.051462747995240005, 0.4934126510932154, 0.8682715826916964, -187.9128456441535 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 118, "378D", //(378_D_N, 378_D_CA, 378_D_C)
[ [ 0.7940017225013606, 0.2988985340513602, 0.5293589812280702, -352.30479211536976 ], [ -0.09355906818972659, 0.9204850078221517, -0.37941277144308616, 18.677767363282182 ], [ -0.6006729271611487, 0.2517280610432642, 0.7588313500765872, -160.70142184159613 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 119, "379M", //(379_M_N, 379_M_CA, 379_M_C)
[ [ -0.41617862171285364, -0.6289826678328361, 0.6566400523841641, -326.1571760286236 ], [ 0.6778042497704947, 0.266779301150814, 0.6851351716782172, -2.495059558096349 ], [ -0.6061161223909544, 0.7302120295115464, 0.3152992834346156, -141.40216020097407 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 120, "380L", //(380_L_N, 380_L_CA, 380_L_C)
[ [ -0.8601112583696242, 0.23506320965134872, -0.45271835692207985, -301.28489411930366 ], [ -0.5092191063992785, -0.4479860368433061, 0.7348499251352947, 26.99689606469308 ], [ -0.030075320509380204, 0.8625855309833536, 0.5050165114473228, -142.61818009679686 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 121, "381L", //(381_L_N, 381_L_CA, 381_L_C)
[ [ 0.48873638287273924, 0.8396962637908673, -0.23674233384025242, -328.2997944356 ], [ -0.856988743538734, 0.4112409403075948, -0.31056590679407675, 52.18722127090341 ], [ -0.16342289161672907, 0.3546703731503148, 0.9205986557154227, -131.53025510960876 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 122, "382A", //(382_A_N, 382_A_CA, 382_A_C)
[ [ 0.5535723264175492, -0.23601374758530025, 0.7986583690010419, -339.39979396205393 ], [ 0.4214979438107122, 0.9065045694585732, -0.0242682696966356, 28.966230976489314 ], [ -0.718259815657262, 0.35006710285525133, 0.6012951527415856, -102.82694686555051 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 123, "383T", //(383_T_N, 383_T_CA, 383_T_C)
[ [ -0.78289085497395, -0.5252443216942665, 0.33346710741254176, -303.0659943211679 ], [ 0.4008584251962819, -0.015939396138216086, 0.9160013420295333, 23.793335215551366 ], [ -0.4758092392412264, 0.8508021733509377, 0.2230274191080694, -91.05003774323178 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 124, "384T", //(384_T_N, 384_T_CA, 384_T_C)
[ [ -0.45542504430681235, 0.6376859115219022, -0.6212444826834429, -296.6611560128718 ], [ -0.8869439937144448, -0.26469192559117716, 0.3785083044540843, 61.85973618694314 ], [ 0.07693101476003722, 0.7233912238708147, 0.6861390210406895, -93.93879693859539 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 125, "385S", //(385_S_N, 385_S_CA, 385_S_C)
[ [ 0.8006612047197028, 0.5026124195354191, 0.32607114405532106, -328.12942156791144 ], [ -0.3880277918690762, 0.8497111812172246, -0.3569668629601882, 66.84321160982218 ], [ -0.4564822756624872, 0.15928485252270763, 0.8753583653354872, -72.22019199197254 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 126, "386R", //(386_R_N, 386_R_CA, 386_R_C)
[ [ -0.10893359315769863, -0.5908351258950956, 0.7994043571874496, -312.0975212378553 ], [ 0.6363196254175654, 0.5763865396540762, 0.5127142393322516, 43.33423624576379 ], [ -0.7636954933677442, 0.564528485476266, 0.3131721291808564, -46.709501149276996 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 127, "387F", //(387_F_N, 387_F_CA, 387_F_C)
[ [ -0.9461137311083102, 0.05810891048675529, -0.3185783456708168, -280.5155605333576 ], [ -0.3086393964053368, -0.45961766204559923, 0.8327624677675421, 65.62393904919867 ], [ -0.09803331471934201, 0.8862138338221648, 0.45278528017964415, -47.09374381749517 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 128, "388R", //(388_R_N, 388_R_CA, 388_R_C)
[ [ 0.3201028146609461, 0.8678700101939879, -0.3799155609500975, -302.7025691939668 ], [ -0.9426909359679553, 0.2519236542494604, -0.21878818906298078, 96.39456046301834 ], [ -0.09416999145163113, 0.42817767087446956, 0.8987746630132127, -38.83634094637204 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 129, "389E", //(389_E_N, 389_E_CA, 389_E_C)
[ [ 0.7137864933018602, -0.08204109324107976, 0.6955415882602866, -322.2411046197347 ], [ 0.2742834269579281, 0.9465335294190148, -0.16983191508605777, 77.58233120180043 ], [ -0.6444202384128298, 0.3119992575398081, 0.6981253609621699, -11.395623695298644 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 130, "390L", //(390_L_N, 390_L_CA, 390_L_C)
[ [ -0.6291763433689096, -0.5574131910354596, 0.5416896375274258, -289.11762673777895 ], [ 0.5418529804823009, 0.18509876013727514, 0.8198376647471685, 65.96984971919062 ], [ -0.557254409123315, 0.8093386086523424, 0.18557623785742433, 3.730264436400352 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 131, "391K", //(391_K_N, 391_K_CA, 391_K_C)
[ [ -0.25260542140724435, -0.07932017515360698, 0.9643126105621919, -274.6869868543816 ], [ -0.9640010921865513, 0.10614621144457746, -0.2437926907417436, 101.95346173548819 ], [ -0.08302045132860265, -0.9911817651720415, -0.10327784395324562, 0.5734210928430423 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 132, "392L", //(392_L_N, 392_L_CA, 392_L_C)
[ [ -0.7261850570601389, 0.14014160193015135, 0.673064331480293, -243.31798699027397 ], [ 0.6607642418873322, -0.1280874031520917, 0.7395838247263204, 85.79002361976775 ], [ 0.1898575243323842, 0.9818115646916187, 0.0004145988104339967, -14.436757278236433 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 133, "393Q", //(393_Q_N, 393_Q_CA, 393_Q_C)
[ [ -0.5206032259035653, -0.512737073978402, 0.6826953743411978, -212.4704158917348 ], [ -0.7951281228398518, 0.5824482190531285, -0.16889446524682958, 106.4814231152374 ], [ -0.31103625101903315, -0.6307572949161235, -0.7109160889036901, -4.921979770514886 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 134, "394H", //(394_H_N, 394_H_CA, 394_H_C)
[ [ -0.43750023446795777, -0.5039015215348679, 0.7447662730248511, -185.43933846930676 ], [ -0.17118998927984963, -0.7664162731414006, -0.6191123353918988, 113.36509590590393 ], [ 0.8827726391381004, -0.3983583221910785, 0.2490444031297672, -31.20293811557237 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 135, "395K", //(395_K_N, 395_K_CA, 395_K_C)
[ [ -0.8967954409632619, 0.33034076132227264, -0.29433470484557184, -158.60769814115736 ], [ 0.2642404633986436, -0.13368927075851164, -0.9551461439941678, 95.43928281837485 ], [ -0.3548730964307927, -0.9343458461903058, 0.032602532670678315, -9.911729156653767 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 136, "396E", //(396_E_N, 396_E_CA, 396_E_C)
[ [ 0.3814469430236681, 0.92113053649973, -0.07756780508444239, -181.44990592704286 ], [ 0.857079933860991, -0.38385817600708044, -0.3436086257435291, 64.51601389415558 ], [ -0.3462834339536587, 0.06458665063491671, -0.9359040270930697, -6.413519796555027 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 137, "397Y", //(397_Y_N, 397_Y_CA, 397_Y_C)
[ [ 0.3999256518336031, -0.20689263322908533, 0.8928913211141671, -188.0535830398132 ], [ 0.12095762319770643, -0.9537554276345478, -0.2751723780614666, 64.41120549367643 ], [ 0.9085310816894651, 0.21805050463866443, -0.3564060199140784, -44.436965359060316 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 138, "398L", //(398_L_N, 398_L_CA, 398_L_C)
[ [ -0.8747613210811525, -0.36822686228195595, 0.3149628692946871, -150.08555877145878 ], [ -0.14652727606694885, -0.4185539998656275, -0.896293649740354, 64.87269321754044 ], [ 0.4618683670795791, -0.8301936684235007, 0.3121795707614024, -51.4667258901668 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 139, "399C", //(399_C_N, 399_C_CA, 399_C_C)
[ [ -0.45946476989030743, 0.7148860917672301, -0.5270958176910743, -144.19273949206348 ], [ 0.6809840795149099, -0.09745573948817193, -0.7257844461601836, 34.71947240749201 ], [ -0.5702217188751035, -0.692416243771468, -0.4420485682416959, -28.599231946387565 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 140, "400V", //(400_V_N, 400_V_CA, 400_V_C)
[ [ 0.7061003990302304, 0.6571834560721503, 0.26368946045379527, -173.87498320619977 ], [ 0.7002752934736535, -0.7033117181168416, -0.1223402652846764, 14.53201068991363 ], [ 0.10505588912449434, 0.27103972444016505, -0.9568180223718935, -43.352608770945544 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 141, "401K", //(401_K_N, 401_K_CA, 401_K_C)
[ [ -0.030306992396889915, -0.4698947265401042, 0.8822020359201488, -161.17858800725372 ], [ -0.16624117424403517, -0.8679450953974908, -0.46801194788317624, 23.131714081067585 ], [ 0.9856192764946408, -0.1608423369179931, -0.051811045703977085, -78.96255247859818 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 142, "402A", //(402_A_N, 402_A_CA, 402_A_C)
[ [ -0.9902739089215081, -0.06362770538927683, -0.12372994955228953, -125.79974588248778 ], [ 0.13552589039880217, -0.23999034295934762, -0.9612686244322481, 11.284856546754026 ], [ 0.03146932380793373, -0.9686878498221964, 0.2462793764526809, -69.39916289958802 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 143, "403M", //(403_M_N, 403_M_CA, 403_M_C)
[ [ 0.14296891599925346, 0.9381276961606143, -0.3153986600389678, -141.2325563127504 ], [ 0.8308285044058391, -0.28693583425513547, -0.47685618721692824, -20.738002314765353 ], [ -0.5378511739549846, -0.19386658483784003, -0.8204461359276168, -54.217418271599804 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 144, "404I", //(404_I_N, 404_I_CA, 404_I_C)
[ [ 0.61855759647308, -0.03757183456951033, 0.7848406571352733, -158.24261816912508 ], [ 0.3251261372374754, -0.897095119699746, -0.29918780238495124, -27.79777891368469 ], [ 0.7153177578744181, 0.44023709913861153, -0.5426894156066946, -88.34125819476917 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 145, "405L", //(405_L_N, 405_L_CA, 405_L_C)
[ [ -0.7007312553204472, -0.5275174049150633, 0.48031353856487047, -123.5733454977798 ], [ -0.06583157335071008, -0.622566763673849, -0.7797928114049527, -27.541662877668607 ], [ 0.7103815254968278, -0.5780449915710256, 0.4015246890945222, -105.47965670072772 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 146, "406L", //(406_L_N, 406_L_CA, 406_L_C)
[ [ -0.8052857882134828, 0.30881336001103293, -0.5061117544377781, -107.68748132946003 ], [ 0.37430195168556973, -0.3972119875231515, -0.8379264203569913, -50.08326233935941 ], [ -0.45979652920151864, -0.86420885533915, 0.20427972510171713, -78.60841811514827 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 147, "407N", //(407_N_N, 407_N_CA, 407_N_C)
[ [ 0.5731981369850326, 0.7913418708966482, -0.2126545064716017, -137.71504272263772 ], [ 0.5646004323677929, -0.5695001864403366, -0.5974076409073817, -74.09921303003364 ], [ -0.5938604613265553, 0.22236812048990864, -0.7732284083392154, -73.41546326347795 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 148, "408S", //(408_S_N, 408_S_CA, 408_S_C)
[ [ 0.12068113366892157, -0.9893398907843185, 0.08150242008220832, -132.60018909748226 ], [ 0.4324053698627525, -0.021514299188345187, -0.9014226151169553, -96.02909848487874 ], [ 0.8935668190405461, 0.14402678720747483, 0.42519951137842654, -104.51445767027977 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 149, "409S", //(409_S_N, 409_S_CA, 409_S_C)
[ [ 0.19173386618291755, 0.6071292773531562, 0.771123962239007, -139.39131047434287 ], [ 0.8799211734785337, 0.24167410130068917, -0.4090627790751249, -132.90701585298856 ], [ -0.43471468001749297, 0.7569594898943659, -0.4878887963831002, -97.55732639769461 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 150, "410M", //(410_M_N, 410_M_CA, 410_M_C)
[ [ -0.90472690275126, -0.3373684586086442, -0.260099508984811, -105.125852908138 ], [ 0.24830320418759882, 0.07848639581460502, -0.9654974906555712, -148.88603312338165 ], [ 0.3461426732263502, -0.9380950958234783, 0.01276091547090059, -102.79459479907119 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 151, "411Y", //(411_Y_N, 411_Y_CA, 411_Y_C)
[ [ 0.2917377747082768, 0.0649888838340869, 0.9542879627168428, -115.40410840550297 ], [ 0.9249575190819825, -0.27326242213326085, -0.26416138352069474, -182.09568227046586 ], [ 0.24360348663754128, 0.9597416806765575, -0.13983292770363545, -117.35097141176331 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 152, "412P", //(412_P_N, 412_P_CA, 412_P_C)
[ [ -0.9746742889250116, 0.2183648381515562, 0.048237205221360324, -81.80985038842269 ], [ 0.15837381848764648, 0.8263016812355344, -0.5405027892665971, -200.95165121972272 ], [ -0.1578852878711787, -0.5191746614063268, -0.8399582768372849, -114.0787694941628 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 153, "413L", //(413_L_N, 413_L_CA, 413_L_C)
[ [ 0.04972167640383509, 0.2738802737121172, 0.9604776679168382, -74.28096297644251 ], [ 0.5032849137771861, -0.8375175490621571, 0.21276430757344555, -208.71799848889663 ], [ 0.862688849157031, 0.4728149222309857, -0.1794825864974767, -150.58411919855843 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 154, "414V", //(414_V_N, 414_V_CA, 414_V_C)
[ [ -0.8397598730222414, -0.011948987963517385, 0.5428264707513043, -36.3000694580669 ], [ 0.0026695060005974625, -0.999836591005484, -0.01787917851147209, -216.0067227476756 ], [ 0.5429514061123533, -0.013565138155582041, 0.8396545466008258, -154.79541450326445 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 155, "415T", //(415_T_N, 415_T_CA, 415_T_C)
[ [ -0.8407773043709664, -0.5411218009414471, 0.01675473069243578, -12.352786488167457 ], [ 0.4928124288896083, -0.7521705879956319, 0.4374646459843699, -231.71608047086238 ], [ -0.22411924144663697, 0.3760672853362682, 0.899079508449166, -128.90497010626027 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 156, "416A", //(416_A_N, 416_A_CA, 416_A_C)
[ [ 0.5584793647011643, -0.24534226042461066, 0.7924064452367398, -17.67755704472826 ], [ -0.07291106088808647, 0.9370451981546112, 0.34151174769773174, -229.49461700518071 ], [ -0.8263079186375727, -0.24850245846898167, 0.5054322424734244, -90.75067852305078 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 157, "417T", //(417_T_N, 417_T_CA, 417_T_C)
[ [ -0.836401882295681, -0.5408945524365865, 0.08868469110656264, 17.44939801457108 ], [ -0.03327340140113831, 0.21160526737529994, 0.9767886626994738, -214.22714393384152 ], [ -0.5471058143096638, 0.8140370347615964, -0.19498444549529337, -83.2601496314482 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 158, "418Q", //(418_Q_N, 418_Q_CA, 418_Q_C)
[ [ -0.3309282786430732, 0.4295053674731828, -0.8402449724372403, 14.088224851319971 ], [ -0.8811248930897034, 0.1780952006130838, 0.43806508910919195, -179.83920181867578 ], [ 0.3377949040053614, 0.8853288873804535, 0.31951113282278126, -100.65641953108077 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 159, "419D", //(419_D_N, 419_D_CA, 419_D_C)
[ [ 0.9679908854622207, 0.25049133938761003, -0.015739585568467527, -23.558444873247065 ], [ -0.23992353331126667, 0.94191917238578, 0.23499993799046998, -174.87515225955661 ], [ 0.07369086663558422, -0.22370150107652023, 0.9718679409213004, -93.1223677270235 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 160, "420A", //(420_A_N, 420_A_CA, 420_A_C)
[ [ 0.18277701741764846, -0.8934205938494986, 0.4103561920935494, -15.958639761774876 ], [ 0.2544915325404048, 0.4461554730599513, 0.8580089473448906, -176.50927013975658 ], [ -0.9496455242716331, -0.05239214010002316, 0.30891494280237536, -55.300561584676785 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 161, "421D", //(421_D_N, 421_D_CA, 421_D_C)
[ [ -0.6539018421540465, -0.4796655585730843, -0.5850925847644246, 5.486696327321974 ], [ -0.6686571438098247, 0.004558419679699718, 0.7435568874283334, -144.4135122634195 ], [ -0.35399153218633894, 0.8774395550286042, -0.32371256140528964, -55.07797993894761 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 162, "422S", //(422_S_N, 422_S_CA, 422_S_C)
[ [ 0.3527130734533169, 0.3833400716615722, -0.8536063948176675, -17.89042344003577 ], [ -0.6951629453678775, 0.7179901816199412, 0.03519344377609317, -124.18378158382184 ], [ 0.6263720677062254, 0.5809823478866327, 0.5197283369622665, -78.0503930753896 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 163, "423S", //(423_S_N, 423_S_CA, 423_S_C)
[ [ 0.9529054412649356, -0.2755517714749534, 0.1266587590523978, -49.00846842704233 ], [ 0.20378009184809925, 0.8910939269341763, 0.40549388102338724, -135.66027342962767 ], [ -0.22459940822285918, -0.3605867920737582, 0.9052802169537881, -58.26613183258546 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 164, "424R", //(424_R_N, 424_R_CA, 424_R_C)
[ [ -0.20210026237531295, -0.9755003931136906, 0.08691643678191475, -32.03635485898899 ], [ -0.0834920791603531, 0.10558598976401337, 0.9908989209213193, -124.9318157632882 ], [ -0.9757994448990587, 0.19300409788548167, -0.10278551227726805, -25.42165270267414 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 165, "425K", //(425_K_N, 425_K_CA, 425_K_C)
[ [ -0.44208131278429863, -0.06447781556939723, -0.8946545278409482, -21.766205099758483 ], [ -0.8834615469609781, 0.20378045977348558, 0.4218639819370958, -90.7796158683397 ], [ 0.1551122429991475, 0.9768910561131706, -0.1470512038640548, -39.83998158978078 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 166, "426L", //(426_L_N, 426_L_CA, 426_L_C)
[ [ 0.7624671023849686, 0.3093470324704903, -0.5682854311720222, -57.23791360957116 ], [ -0.35405564325439937, 0.934616115518231, 0.033724147036302235, -85.05843424543976 ], [ 0.541561186995877, 0.17549111121454, 0.8221401040116597, -53.49488267016712 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 167, "427A", //(427_A_N, 427_A_CA, 427_A_C)
[ [ 0.6625119345115749, -0.6691342458014494, 0.3366560525602595, -72.46364389270116 ], [ 0.29372801914143887, 0.6455220061000054, 0.7050001350438673, -94.81207345717722 ], [ -0.6890586240669058, -0.3681856878483329, 0.6242095095900676, -19.4885072133021 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 168, "428H", //(428_H_N, 428_H_CA, 428_H_C)
[ [ -0.5392036321053965, -0.7827462139955069, -0.31075361236846205, -48.51880195843572 ], [ -0.46165524789570045, -0.033892789867279415, 0.8864117050701508, -70.14212903133884 ], [ -0.7043677130694647, 0.6214174468670095, -0.3430837820628947, -1.7887467554970158 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 169, "429L", //(429_L_N, 429_L_CA, 429_L_C)
[ [ 0.021567913246255364, 0.2680689497583127, -0.9631582753076865, -59.29601654028437 ], [ -0.849076228166107, 0.513531694151036, 0.12391431662880556, -42.13548229016554 ], [ 0.5278298815730491, 0.8151222222942042, 0.23868677977788985, -26.27962049092804 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 170, "430L", //(430_L_N, 430_L_CA, 430_L_C)
[ [ 0.9847415183653015, -0.005754495013739748, -0.173928225986288, -96.51699663050678 ], [ 0.03460056210591747, 0.9859730001274963, 0.16327903760599083, -50.366963432063756 ], [ 0.17054894637480164, -0.1668056617945592, 0.9711276579748536, -21.082850816644928 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 171, "431N", //(431_N_N, 431_N_CA, 431_N_C)
[ [ 0.1292829913887037, -0.9689235681900985, 0.21088628960498068, -91.83915857728 ], [ 0.16676574427641444, 0.23088702142599793, 0.9585824794314594, -50.98130583082977 ], [ -0.9774840636417158, -0.08875980139003352, 0.19143302479896174, 17.144724514780016 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 172, "432A", //(432_A_N, 432_A_CA, 432_A_C)
[ [ -0.5685169033001597, -0.4309159663671647, -0.7007852457007422, -74.1984027531644 ], [ -0.7670636887914781, -0.030205912495388225, 0.6408595011294843, -16.669217334951753 ], [ -0.2973244490444865, 0.9018863746504732, -0.31336741888750813, 16.655789839082072 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 173, "433V", //(433_V_N, 433_V_CA, 433_V_C)
[ [ 0.5143903400851967, 0.295444648024053, -0.8050559222693658, -102.91511052456262 ], [ -0.6383374731838055, 0.7588031184757981, -0.1293951224765381, -1.3954469930028184 ], [ 0.5726498479492944, 0.5804569642491573, 0.5789143842558264, -4.2975508142028165 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 174, "434T", //(434_T_N, 434_T_CA, 434_T_C)
[ [ 0.8464562526238653, -0.5084187865065862, 0.15818391170132748, -128.04990658149825 ], [ 0.3830437146795509, 0.7877939096624134, 0.4823466269637428, -18.012990229101355 ], [ -0.3698504090013339, -0.34769396518487267, 0.8615797012091027, 19.93637035372771 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 175, "435D", //(435_D_N, 435_D_CA, 435_D_C)
[ [ -0.320205186669698, -0.9445757112092004, -0.07242488676868725, -109.3049904218699 ], [ -0.15642592671847966, -0.02268325740266236, 0.9874291869718402, -3.107279293887295 ], [ -0.9343444589015814, 0.32750907716769906, -0.1404928343107699, 50.20786801826097 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 176, "436A", //(436_A_N, 436_A_CA, 436_A_C)
[ [ -0.2008332285815266, 0.08250725547425304, -0.9761447469979165, -109.90000628602441 ], [ -0.9256760810113472, 0.31013045977867304, 0.21666308167508028, 30.704144973017335 ], [ 0.32060849542860514, 0.9471069901080739, 0.014090491383584392, 31.772869139790245 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 177, "437L", //(437_L_N, 437_L_CA, 437_L_C)
[ [ 0.9337808025987168, 0.12951630034933387, -0.3335849826384837, -148.0152158464157 ], [ -0.13124694548066973, 0.9911961604478268, 0.017447372738613184, 27.109281685258313 ], [ 0.332907873142241, 0.027489988310422046, 0.9425585650465452, 27.373571993338825 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 178, "438V", //(438_V_N, 438_V_CA, 438_V_C)
[ [ 0.3716223992626223, -0.8781161631758165, 0.30134497894551837, -152.2078961830938 ], [ 0.3500453617127548, 0.4331666986870258, 0.8305629752715745, 17.626308772071635 ], [ -0.8598633828171071, -0.2031713934538181, 0.468354937811437, 64.61481506362006 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 179, "439W", //(439_W_N, 439_W_CA, 439_W_C)
[ [ -0.6083005658033606, -0.5607791333199945, -0.5616913612262406, -130.55775592146009 ], [ -0.6117441817130397, -0.11964314384293033, 0.781955608888141, 48.07522099690131 ], [ -0.5057069089735228, 0.8192754614683976, -0.27027437994044023, 73.9816022506306 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 180, "440V", //(440_V_N, 440_V_CA, 440_V_C)
[ [ 0.32906586075727673, 0.35899064051369683, -0.8734079111776113, -152.63974977178628 ], [ -0.7796396768081865, 0.6251468771177617, -0.03678799228466607, 71.61702237526677 ], [ 0.5328016832091459, 0.6930491339389119, 0.48559784216574226, 52.343470082109334 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 181, "441I", //(441_I_N, 441_I_CA, 441_I_C)
[ [ 0.9512954975254706, -0.2932077062941039, 0.0952161612202362, -184.25858725847337 ], [ 0.25898043039595176, 0.9276369410040614, 0.2691074921969945, 56.88951015572633 ], [ -0.16723041906212865, -0.23134162326390442, 0.9583919032869208, 68.62852200195528 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 182, "442A", //(442_A_N, 442_A_CA, 442_A_C)
[ [ -0.1566621119223425, -0.9819866185174452, 0.1056374163859803, -169.13094466257715 ], [ 0.13267228047983998, 0.08506484411485915, 0.9875029307743818, 60.99826381630842 ], [ -0.9787006941247717, 0.16871945160054383, 0.11695596595001743, 103.89149656339131 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 183, "443K", //(443_K_N, 443_K_CA, 443_K_C)
[ [ -0.4826253861279263, -0.11588507655826936, -0.8681263650505923, -160.1986979344707 ], [ -0.8748850775083139, 0.017834656576967915, 0.48400209315461995, 97.8570376408071 ], [ -0.04060588403349695, 0.9931024992697711, -0.10999358220365646, 95.80914973561156 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 184, "444S", //(444_S_N, 444_S_CA, 444_S_C)
[ [ 0.7990933083939522, 0.3986803121648305, -0.4500043257260548, -197.80983573045924 ], [ -0.3909944977076225, 0.9132116534742145, 0.1147509421366777, 104.40935558845611 ], [ 0.45669813579909857, 0.08425250531019539, 0.8856231298394276, 89.58237099811736 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 185, "445G", //(445_G_N, 445_G_CA, 445_G_C)
[ [ 0.008954802816788173, 0.6286298228836551, -0.7776531085823402, -203.57774815275775 ], [ -0.45359021842830793, -0.6905195201907623, -0.5634169912078258, 104.44049527311387 ], [ -0.8911653748057498, 0.3577811314432105, 0.2789568725278942, 127.93705271699034 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 186, "446I", //(446_I_N, 446_I_CA, 446_I_C)
[ [ 0.9233468065165064, -0.2658422056806158, 0.2770534182691751, -230.38982593506253 ], [ 0.38243845376871743, 0.5724202900175922, -0.7253108579465554, 76.61496759948831 ], [ 0.03422724024460553, 0.7756693453109004, 0.6302107288599562, 127.17873147434443 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 187, "447S", //(447_S_N, 447_S_CA, 447_S_C)
[ [ -0.05443008334236132, 0.9397206369672072, 0.33758330895541916, -233.02268955881587 ], [ 0.43072547591890453, 0.3271089027461647, -0.8411155272247881, 49.261051469582966 ], [ -0.9008401247845129, 0.09962374316470211, -0.42256618342701857, 154.2525115933399 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 188, "448S", //(448_S_N, 448_S_CA, 448_S_C)
[ [ -0.7325358513799477, 0.1904267283091027, -0.653550983159345, -214.2237039413758 ], [ 0.6552804829172985, 0.4573123512236923, -0.6012261655366726, 15.702793136464152 ], [ 0.18438740507618023, -0.8686789248990386, -0.4597806110369631, 149.88723305346622 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 189, "449Q", //(449_Q_N, 449_Q_CA, 449_Q_C)
[ [ 0.6757184295303835, 0.6837700550991496, -0.27543259745843124, -247.13440699645997 ], [ 0.7282266442489052, -0.5611906935794865, 0.3933839854973173, -2.8958001932581743 ], [ 0.11441397903652717, -0.4663941650465633, -0.8771464667953398, 141.84749969257015 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 190, "450Q", //(450_Q_N, 450_Q_CA, 450_Q_C)
[ [ 0.5193603669480634, -0.3807944256965026, 0.7650231464485526, -255.47414800684632 ], [ -0.5446816887653503, -0.8373243614369632, -0.04700820851625625, 23.593238517351445 ], [ 0.6584729813494999, -0.39227989892760556, -0.642284838471257, 115.1349765580559 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 191, "451Q", //(451_Q_N, 451_Q_CA, 451_Q_C)
[ [ -0.8682819280318704, -0.493652411242587, -0.04892637660442082, -219.49760438957466 ], [ -0.18547714337262833, 0.4145355812227209, -0.8909312438042963, 20.05665124540502 ], [ 0.46009208071818336, -0.7645047735460654, -0.4514949927581611, 101.17581426750463 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 192, "452S", //(452_S_N, 452_S_CA, 452_S_C)
[ [ -0.16612521537251282, 0.6235970174366939, -0.7638908120022755, -225.5157830082888 ], [ 0.9712557688313366, 0.23739991492422963, -0.017421593090698385, -17.66140630020307 ], [ 0.17048356029038564, -0.7448275238187886, -0.6451103126076084, 95.65745188021158 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 193, "453M", //(453_M_N, 453_M_CA, 453_M_C)
[ [ 0.9327591723562313, 0.28017052070684506, 0.22685855882503878, -259.907649910668 ], [ 0.11820073507768612, -0.8321858007772448, 0.5417521381700587, -10.089551162223133 ], [ 0.34057145008413864, -0.47850942760993864, -0.8093453002742377, 80.0161353665939 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 194, "454R", //(454_R_N, 454_R_CA, 454_R_C)
[ [ -0.12386349450598759, -0.8086806234780989, 0.5750595481685687, -243.6880446715585 ], [ -0.8098789023208303, -0.2524734681711348, -0.529484004900286, 17.182096135877774 ], [ 0.5733707337354884, -0.5313123347718479, -0.6236611296328444, 57.99907245310514 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 195, "455L", //(455_L_N, 455_L_CA, 455_L_C)
[ [ -0.8501968922111575, -0.06919201750682594, -0.5218981789466453, -215.94011845801023 ], [ 0.3703960154682264, 0.6258637806806326, -0.6863682100428358, -6.362413614817703 ], [ 0.37412836861129745, -0.7768571250529147, -0.5064789936954632, 46.041451053633054 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 196, "456A", //(456_A_N, 456_A_CA, 456_A_C)
[ [ 0.4046637300411627, 0.7385478098570728, -0.5392535564504827, -241.7314831768892 ], [ 0.897126501127935, -0.20633275281999003, 0.39062877017402853, -33.84274282530285 ], [ 0.17723235190879266, -0.6418519515192315, -0.7460661939585737, 38.249114653815894 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 197, "457N", //(457_N_N, 457_N_CA, 457_N_C)
[ [ 0.8013292813371949, -0.18885513888484953, 0.5676311473028932, -263.9731605341555 ], [ -0.3945140721889682, -0.88011378203307, 0.26411811282137226, -9.351213180929193 ], [ 0.4496999329736354, -0.4355840529590469, -0.7797669543467959, 18.27549738569708 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 198, "458L", //(458_L_N, 458_L_CA, 458_L_C)
[ [ -0.5448664665479849, -0.7642146752284634, 0.3451035551786852, -233.52856064757435 ], [ -0.6041160990403134, 0.07234412620791106, -0.7936057373050728, 2.867819149855656 ], [ 0.5815189356433959, -0.6408917674371142, -0.5010922768507853, -2.0268391824434717 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 199, "459L", //(459_L_N, 459_L_CA, 459_L_C)
[ [ -0.6143631709394437, 0.305740980304044, -0.7273790945277119, -220.7681458021381 ], [ 0.747981715083741, 0.5191264736186084, -0.4135590142757058, -32.727774921168916 ], [ 0.25115980588783904, -0.7981416900219098, -0.5476208492701429, -9.25182585873266 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 200, "460M", //(460_M_N, 460_M_CA, 460_M_C)
[ [ 0.8134874554905419, 0.5765553716189256, -0.07630244568089604, -256.1060724529384 ], [ 0.5327321333527506, -0.6860846955274835, 0.49546368651636546, -41.82693316610657 ], [ 0.233312309690146, -0.4437022582997778, -0.865270866334194, -22.034051495968953 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 201, "461L", //(461_L_N, 461_L_CA, 461_L_C)
[ [ 0.18925954122374863, -0.8285998201389251, 0.5268805975944807, -251.27404080309924 ], [ -0.717601186464925, -0.482973720927474, -0.5017817474536138, -15.043848694398829 ], [ 0.6702457483937575, -0.28312315864162624, -0.6860115988821767, -49.4535039373969 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 202, "462L", //(462_L_N, 462_L_CA, 462_L_C)
[ [ -0.8235239118092995, -0.14900229528402908, -0.5473633917959804, -224.47125734953204 ], [ 0.3307385519459355, 0.6578165194949042, -0.6766752817536015, -38.22465128803849 ], [ 0.46089085143341113, -0.7382924506452042, -0.4924468300083614, -64.85371760543347 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 203, "463S", //(463_S_N, 463_S_CA, 463_S_C)
[ [ 0.37973842477070596, 0.6897245586246572, -0.6165052813907476, -251.56601261933903 ], [ 0.8988393753501484, -0.11744546371008492, 0.4222491449062726, -64.57950316498105 ], [ 0.218829856447465, -0.7144834471728124, -0.6645501468235763, -72.78527577862262 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 204, "464H", //(464_H_N, 464_H_CA, 464_H_C)
[ [ 0.8235363215977547, -0.2043053314873368, 0.5291949154471243, -274.4113864801385 ], [ -0.41460551409484886, -0.8534808767712249, 0.3157097728420321, -38.65925861971468 ], [ 0.3871565506244026, -0.47940559499412483, -0.7875786188038268, -89.70580787483563 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 205, "465V", //(465_V_N, 465_V_CA, 465_V_C)
[ [ -0.5895635776348142, -0.7708439654565089, 0.24127654018930839, -244.9131653693022 ], [ -0.5901698534802472, 0.20715927582772412, -0.7802464857220676, -24.21730693716211 ], [ 0.551465621747668, -0.6023990499310382, -0.5770624339468309, -109.98837461148722 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 206, "466R", //(466_R_N, 466_R_CA, 466_R_C)
[ [ -0.5340142920527218, 0.35437200085155035, -0.7676257036446208, -234.5111219170352 ], [ 0.7871556943781922, 0.5397466154610311, -0.29842838991347304, -60.04824097253842 ], [ 0.3085687098385324, -0.7636059691430067, -0.5671783451417666, -120.30210813893143 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 207, "467H", //(467_H_N, 467_H_CA, 467_H_C)
[ [ 0.7972288806535474, 0.5715014072542796, -0.19445373063603066, -270.3269027838067 ], [ 0.5583578445831868, -0.5756274186294729, 0.5974023705296017, -68.43712946470356 ], [ 0.22948339644582075, -0.5848411890661838, -0.7780090965582231, -131.06144871787802 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 208, "468A", //(468_A_N, 468_A_CA, 468_A_C)
[ [ 0.4038329452495043, -0.5578596575453899, 0.7250596905182712, -274.58269182699854 ], [ -0.7521099074209697, -0.6536592632881992, -0.08402532164025904, -36.23444661365524 ], [ 0.5208163202995191, -0.5113923835968512, -0.6835409208736517, -151.81602311585522 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 209, "469S", //(469_S_N, 469_S_CA, 469_S_C)
[ [ -0.897865740310366, -0.4290515832495451, -0.09875146220677991, -239.389807675022 ], [ -0.17898922710256854, 0.5606475327497986, -0.808478324139113, -41.64424329272799 ], [ 0.4022436686364915, -0.7082295411324088, -0.5801818233185475, -166.60166120427553 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 210, "470N", //(470_N_N, 470_N_CA, 470_N_C)
[ [ -0.086555802262362, 0.7268948374004097, -0.6812723306104179, -248.0496687473988 ], [ 0.9482958623675707, 0.2697160315554211, 0.16729680133983588, -77.33488355337846 ], [ 0.3053572506282834, -0.6315672234074745, -0.7126568541774567, -177.9989253380022 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 211, "471K", //(471_K_N, 471_K_CA, 471_K_C)
[ [ 0.9049621538864434, 0.2496552336270212, 0.34455154092798146, -280.13741594548685 ], [ -0.0266010528096866, -0.7749917916372755, 0.6314112026914489, -64.34800657639546 ], [ 0.42465972733781254, -0.5805686757115675, -0.6946971489504556, -195.24875061521925 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 212, "472G", //(472_G_N, 472_G_CA, 472_G_C)
[ [ -0.2352264169202011, -0.7381872069616502, 0.6322564197072589, -259.9311656517801 ], [ -0.8761723232253497, -0.12052090034757419, -0.4666870178109775, -35.69384939999705 ], [ 0.4207024991567999, -0.6637426911519011, -0.6184295005460645, -211.19282856307603 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 213, "473M", //(473_M_N, 473_M_CA, 473_M_C)
[ [ -0.8594566981053804, 0.10772374247674564, -0.49972970632993274, -233.8868451748582 ], [ 0.4048101204920628, 0.740392221449599, -0.5366079805259418, -61.21246608370721 ], [ 0.31219056748883156, -0.6634869657527055, -0.6799427151218472, -223.53816933892034 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 214, "474E", //(474_E_N, 474_E_CA, 474_E_C)
[ [ 0.48793953023525827, 0.7743682296400682, -0.40282609120798246, -261.09895414649355 ], [ 0.7981104581746945, -0.20891911515460154, 0.565130515788342, -84.11225955234559 ], [ 0.3534610464902386, -0.597249234614097, -0.7199712774596642, -238.44966402386416 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 215, "475H", //(475_H_N, 475_H_CA, 475_H_C)
[ [ 0.6896083568211497, -0.22514776408847334, 0.68829412210798, -277.92299282871926 ], [ -0.5665681080121185, -0.7597024564530951, 0.31914378678379646, -54.203216889927695 ], [ 0.4510442253105587, -0.6100497209122336, -0.6514587053903752, -255.70547977126694 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 216, "476L", //(476_L_N, 476_L_CA, 476_L_C)
[ [ -0.6566328732805429, -0.6789512699580682, 0.32841809138606365, -244.61476867565008 ], [ -0.6309695061921137, 0.25597464069383963, -0.7323622502405076, -40.9556680526361 ], [ 0.41317157693023704, -0.6881149296040496, -0.5964789113397012, -270.07122502059786 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 217, "477L", //(477_L_N, 477_L_CA, 477_L_C)
[ [ -0.4814659288290961, 0.5363221278108912, -0.6932165135057977, -235.84217735825266 ], [ 0.816368178896052, 0.5622426767663471, -0.13200821530706597, -76.48916074903872 ], [ 0.3189569812101305, -0.6294773607072863, -0.7085370113792994, -282.4172503896465 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 218, "478N", //(478_N_N, 478_N_CA, 478_N_C)
[ [ 0.8200001959222072, 0.5717390605052162, 0.026723124445964867, -270.2715536797652 ], [ 0.3838126646430813, -0.583908153559723, 0.7153594248110893, -78.7681889023462 ], [ 0.4246027757176341, -0.5763381948986962, -0.6982454926125861, -299.3767346223341 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 219, "479M", //(479_M_N, 479_M_CA, 479_M_C)
[ [ 0.23464019467795566, -0.5734607765015183, 0.7849119166223474, -266.0268034621846 ], [ -0.8659703728677226, -0.4901507149515572, -0.09923502379611147, -44.006678526284006 ], [ 0.4416325309087053, -0.6564259398033817, -0.6116091016298061, -315.78375683937935 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 220, "480K", //(480_K_N, 480_K_CA, 480_K_C)
[ [ -0.9387938660338305, -0.2480710860849945, -0.23901216150201515, -231.11061696928505 ], [ 0.0032161910257858207, 0.6874936404335059, -0.7261832760941077, -54.60947094965925 ], [ 0.3444644150163262, -0.6825051139823763, -0.6446169685754197, -328.20963321860205 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 221, "481C", //(481_C_N, 481_C_CA, 481_C_C)
[ [ 0.023820173320541736, 0.76167886661962, -0.6475167206242866, -245.4144125531536 ], [ 0.9407747340980935, 0.2020324274776048, 0.27226053318497656, -87.61475792204077 ], [ 0.33819446924164903, -0.6156526637581919, -0.711755785773324, -342.3133091786479 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 222, "482K", //(482_K_N, 482_K_CA, 482_K_C)
[ [ 0.8667290814503914, 0.17437995173419077, 0.46730325464449957, -272.94894360704063 ], [ -0.19049217218467185, -0.7501458412130985, 0.6332408303695034, -67.32365864895007 ], [ 0.46097009849286974, -0.6378658552892404, -0.6169552001171558, -360.32455069571336 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 223, "483N", //(483_N_N, 483_N_CA, 483_N_C)
[ [ -0.7703563600420699, 0.5109399593129026, -0.3814333972269107, -244.8186971107042 ], [ -0.3180932612795561, 0.2104952510119017, 0.9243962496840695, -44.761426307607465 ], [ 0.5526009008961792, 0.8334459234280086, 0.0003702026002145075, -375.28997781858243 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 224, "484V", //(484_V_N, 484_V_CA, 484_V_C)
[ [ 0.7631982712769133, 0.03385586501631286, 0.645276823637676, -264.8734625775745 ], [ -0.646064011180468, 0.05758183307293868, 0.76110815654371, -13.431613380056191 ], [ -0.011388247333782194, -0.997766562330066, 0.06581941141264615, -364.19079492903165 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 225, "485V", //(485_V_N, 485_V_CA, 485_V_C)
[ [ -0.369074752162349, 0.07172449072277128, 0.9266279861664333, -233.66637656131974 ], [ -0.9020513134020374, 0.21243995873546417, -0.37572954624587834, 6.5161795200903185 ], [ -0.22380182149830694, -0.9745382811172855, -0.01370705406010301, -353.43617171406817 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 226, "486P", //(486_P_N, 486_P_CA, 486_P_C)
[ [ -0.9140472320275376, -0.26970409997894057, 0.30294777780559207, -195.49440549402132 ], [ 0.28643183415137285, -0.9580338571068299, 0.011310748060616145, 4.202050073817921 ], [ 0.2871836729472572, 0.09711244560590533, 0.9529400353122766, -358.04398269252414 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 227, "487V", //(487_V_N, 487_V_CA, 487_V_C)
[ [ -0.15842290139781234, -0.1436029079749747, 0.9768727599507676, -178.11045254662577 ], [ 0.13481131004294547, 0.9769559631139035, 0.1654779647587256, -7.749540105976053 ], [ -0.9781247849825726, 0.157908995808713, -0.13541290206431586, -326.10512278546645 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 228, "488Y", //(488_Y_N, 488_Y_CA, 488_Y_C)
[ [ -0.8192390780222194, 0.08173714192878213, 0.5675970160867815, -143.61778623754108 ], [ -0.5502662252002283, -0.3906525859514597, -0.7379685890966695, 8.750151504104679 ], [ 0.16141379880665585, -0.9169021740179001, 0.3650150528896759, -322.06421003470706 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 229, "489D", //(489_D_N, 489_D_CA, 489_D_C)
[ [ -0.6087925593856026, -0.788658349696286, 0.08596294021878356, -113.68369749013462 ], [ 0.7639224653809694, -0.6120029221353201, -0.2046335509736634, -15.294630855671281 ], [ 0.21399552921261958, -0.05891036201009068, 0.9750566561615027, -318.2605546072089 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 230, "490L", //(490_L_N, 490_L_CA, 490_L_C)
[ [ -0.42569564692702094, 0.28149420632364214, -0.8599675738035755, -108.89978547142553 ], [ 0.00971432892385176, -0.9489042197754608, -0.31541466913554156, -10.173760635800166 ], [ -0.9048142616034213, -0.14262465950370132, 0.40120986839870837, -280.27832306640545 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 231, "491L", //(491_L_N, 491_L_CA, 491_L_C)
[ [ 0.8855611577974899, 0.41086365894649746, -0.21673137649511084, -146.5525164612822 ], [ -0.04742499640138019, -0.38416309207024296, -0.9220464133693917, -12.843157057374272 ], [ -0.4620955588584447, 0.8268069741167683, -0.3207147050497925, -272.8335758942213 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 232, "492L", //(492_L_N, 492_L_CA, 492_L_C)
[ [ 0.358652346357445, -0.7691325019858124, 0.5289647330790056, -147.52304616567127 ], [ 0.6808325773192985, -0.17214185262709777, -0.7119228780105834, -43.76507165279063 ], [ 0.6386199935118386, 0.615469233157235, 0.461911384276015, -296.0060688217484 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 233, "493E", //(493_E_N, 493_E_CA, 493_E_C)
[ [ -0.7748927841297837, -0.4780009388167122, -0.4135895013100965, -118.32387652061973 ], [ 0.5820617579368996, -0.7947389905000366, -0.17202920370212965, -61.003787442758586 ], [ -0.24646558187909023, -0.3740388208051733, 0.8940635757482109, -277.42489497488657 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 234, "494M", //(494_M_N, 494_M_CA, 494_M_C)
[ [ 0.08752788075784243, 0.4715288174076734, -0.8774961221818374, -134.46862259341745 ], [ -0.14950472222128433, -0.8646818476548881, -0.4795556697294545, -53.62641973153176 ], [ -0.9848792860668649, 0.1731643054738767, -0.005187984789143005, -243.19159558303633 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 235, "495L", //(495_L_N, 495_L_CA, 495_L_C)
[ [ 0.9899359697496916, -0.04804035514653824, 0.13311235882941955, -168.47778932252768 ], [ 0.12006946801081836, -0.2127184996174129, -0.9697082874618097, -66.55838283493404 ], [ 0.07490059176893554, 0.9759318440331328, -0.20480951431696742, -256.1544141983162 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 236, "496N", //(496_N_N, 496_N_CA, 496_N_C)
[ [ -0.10344013432512175, -0.9349712795776892, 0.3393064175279615, -153.23325120392965 ], [ 0.8028393816762964, -0.27986486824392753, -0.5264262367629782, -98.73510711429881 ], [ 0.5871533580253313, 0.21795495380391705, 0.779581023545282, -271.2582877563882 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 237, "497A", //(497_A_N, 497_A_CA, 497_A_C)
[ [ -0.6915766616857815, 0.10605571040362254, -0.7144745672886705, -135.58518186276737 ], [ 0.2924574520233316, -0.8633379599835116, -0.411237407840689, -106.73344241467706 ], [ -0.6604470908161916, -0.49335560515949495, 0.566047601434842, -237.86293042994092 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 238, "498H", //(498_H_N, 498_H_CA, 498_H_C)
[ [ 0.6871844447912003, 0.6195721348037732, -0.3793519587555925, -170.10319787719615 ], [ 0.09088835203175315, -0.5913933519794549, -0.8012447882510416, -113.01587313194722 ], [ -0.7207751704255368, 0.5161242805848215, -0.46270820253031886, -221.60633910574595 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 239, "499V", //(499_V_N, 499_V_CA, 499_V_C)
[ [ 0.5896216330206722, -0.581686125084732, 0.5603459482836878, -182.26582465098636 ], [ 0.5124667701531748, -0.26680459932322215, -0.8162065395895503, -136.5891482542689 ], [ 0.6242788954968836, 0.7684117110403019, 0.14078104621470386, -249.95810461519895 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 240, "500L", //(500_L_N, 500_L_CA, 500_L_C)
[ [ -0.7448250068885562, -0.665938440256833, -0.04197264468387374, -150.6618336013127 ], [ 0.6284891968600576, -0.6790287139400144, -0.3793696549213839, -158.068136507276 ], [ 0.22413620533876738, -0.30894335958758934, 0.9242926820131551, -245.77996907750443 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 241, "501R", //(501_R_N, 501_R_CA, 501_R_C)
[ [ -0.18059639951223613, 0.44483914407369096, -0.8772132445323673, -158.6270674886303 ], [ 0.14962799242594524, -0.8690713616867717, -0.47151503918590093, -160.87447462817838 ], [ -0.9721092553647446, -0.21640957510167258, 0.09039077075960511, -207.98690575651472 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ],
     [ 242, "502G", //(502_G_N, 502_G_CA, 502_G_C)
[ [ 0.9821516408723061, 0.18492312105818748, 0.03437431642863064, -195.7408587545467 ], [ 0.10924845870890823, -0.4120877190054338, -0.9045708850690624, -169.92842234418043 ], [ -0.1531108376358926, 0.8921811201448939, -0.4249351953597604, -213.63349380019105 ], [ 0.0, 0.0, 0.0, 1.0 ] ] ]
   ]
 ]

];
