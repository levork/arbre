// This is a partial implementation of Jason Weber and Joseph Penn's
// paper, "Creation and Rendering of Realistic Trees", from SIGGRAPH
// 1995. This code should not be treated as an example on how to
// implement that paper properly (even though there are some confusing
// sections in the paper that I've tried to interpret in a useful
// fashion), but more as a basic example on how one could set up
// nested object instance hierarchies to maximize memory efficiency in
// PRMan 17.0 (as described in my User's Group Meeting talk at
// SIGGRAPH 2012).
// 
// Sections 4.1 through 4.5 of the Penn/Weber paper are mostly
// implemented, 4.6 through 4.9 are not. I also did not implement
// ternary stem splits, handling of -nCurveV (helical stems), or trunk
// lobes.
//
// To use this plugin, compile it as a DSO. On Linux, this can be done
// using a command resembling:
//
// g++ -I$RMANTREE/include -fPIC -o arbre.so -shared arbre.cpp
//
// The resulting shared object can be loaded with RiProcedural2
// using a recent version of prman:
//
// Procedural2 "DynamicLoad" "SimpleBound" "float[6] bound" [-20 20
//     -20 20 0 20] "string dsoname" ["arbre"] "int Seed" [66] "int
//     MaxLevel" [5] "int StartLevel" [0] "int InstanceLevel" [4] "int
//     Leaves" [1] "int InstanceLeaves" [1]
//
// The useful parameters are:
//
// "string Type": tree type. Currently supported: oak, tupelo, and aspen.
//    The parameters are straight out of the paper, slightly modified
//    in a few cases
//
// "int Seed": random seed for the tree. If you are growing a forest
//    of trees, every tree should have its own seed
//
// "int MaxLevel": maximum level of the tree - this overrides the Levels
//    parameter in the paper
//
// "int InstanceLevel": level at which to use instanced branches. For
//    example, if InstanceLevel == MaxLevel, the terminal level of the
//    tree will use instanced branches. For two levels of terminal instanced
//    branches, set InstanceLevel == MaxLevel - 1, and so on. To disable
//    instanced branches, set this parameter to 0.
//
// "int Leaves": toggle switch for whether leaves will be grown or not
//
// "int InstanceLeaves": whether leaves will be instanced. If 1,
//    the procedural requires a "leaf" ObjectInstance to have
//    been previously defined. If 0, the procedural requires a "leaf"
//    inline archive to be previously defined.
//
// "int StartLevel": by default, the procedural starts at Level 0, this
//    parameter overrides that default. This is mainly useful for
//    visualization purposes
//
//
// There are several known bugs with this plugin that ought to be
// addressed before even considering using this for production work.
//
// - If a stem has high curvature, its children may be detached at the
// base. This is "fixed" by switching to polygons. The proper fix for
// subdivision meshes is to interpolate using the cubic b-spline
// rather than just perform linear interpolation as is done here
//
// - The actual formulae for the leaf/branch frequency/stem
// length/stem radius are dependent not only parent properties but
// also grandparent properties. So, the simple-minded instancing being
// used here where we just look for the branch of the right length
// isn't close to correct as far the Penn/Weber model is
// concerned. I've attempted to compensate for this by associating an
// arbitrary parent offset with each stem length, and this seems to be
// "good enough" at a first glance. A more rigorous derivation,
// coupled perhaps by searching for appropriate instances in multiple
// dimensions (e.g. radius and length), may provide more accurate
// results.
//
// - The use of erand48 throughout is completely ad hoc and wacky.
//
// Naturally, there are other useful features that should be implemented
// as well (like switching amongst multiple leaf types).
// 
//
// Julian Fong <jfong@pixar.com>
#include <vector>
#include <deque>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define DLL /* Must be done on Windows to get the correct declspec */
#include "ri.h"
#include "RixInterfaces.h"

#ifndef M_PI
#define M_PI 3.1415926535
#endif

extern "C" {
    PRMANEXPORT RtVoid Subdivide2(RtContextHandle ctx, RtFloat detail,
                                  RtInt n, RtToken toks[], RtPointer vals[]);
}

struct arbreParameters {
    const struct treeParameters *tree;
    float maxStemLengths[5];
    int maxLevel;
    int startLevel;
    int instanceLevel;
    int leaves;
    int instanceLeaves;
    float terminalScale;
};

static void growstem(const struct arbreParameters &tree,
    float baseradius, float length, int level, float parentV,
    unsigned short xsubi[3]);

static void
rotateMatrix(float angle, float axisx, float axisy, float axisz, RtMatrix m)
{
    RtMatrix tm;
    float a, b, c, d;
    float cost, sint, radians, aa, bb, cc, ab, ac, bc, cost1;
    float asint, bsint, csint, abcost1, accost1, bccost1;
    d = sqrtf(axisx * axisx + axisy * axisy + axisz * axisz);
    d = 1.0f/d;
    a = axisx * d;
    b = axisy * d;
    c = axisz * d;
    radians = angle * 0.0174532925f;
    cost = cosf(radians);
    sint = sinf(radians);
    cost1 = 1.0 - cost;
    asint = a*sint; bsint = b*sint; csint = c*sint;
    aa = a*a; bb = b*b; cc = c*c;
    ab = a*b; ac = a*c; bc = b*c;
    abcost1 = ab*cost1; accost1 = ac*cost1; bccost1 = bc*cost1;
    tm[0][3] = tm[1][3] = tm[2][3] = tm[3][0] = tm[3][1] = tm[3][2] = 0.;
    tm[3][3] = 1.0;
    tm[0][0] = aa+(1-aa)*cost;
    tm[0][1] = abcost1+csint;
    tm[0][2] = accost1-bsint;
    tm[1][0] = abcost1-csint;
    tm[1][1] = bb+(1-bb)*cost;
    tm[1][2] = bccost1+asint;
    tm[2][0] = accost1+bsint;
    tm[2][1] = bccost1-asint;
    tm[2][2] = cc+(1-cc)*cost;

    RtMatrix t;
    float *ap, *bp;
    int i, j;
    float *rp = t[0];
    for(i=0;i!=4;i++)for(j=0;j!=4;j++){
        float sum;
        int k;
        ap = &tm[i][0];
        bp = &m[0][j];
        sum=0.;
        for(k=0;k!=4;k++){
            sum += *ap++ * *bp;
            bp +=4;
        }
        *rp++ = sum;
    }
    memcpy(m, t, sizeof(RtMatrix));
}

struct stemParameters {
    float downAngle, downAngleV;
    float rotate, rotateV, branches;
    float length, lengthV, taper;
    float segSplits, splitAngle, splitAngleV;
    int curveRes;
    float curve, curveBack, curveV;
};

struct leafParameters {
    float leaves;
    float leafScale, leafScaleX;
    float attractionUp;
};

struct treeParameters {
    enum Shape {
        k_Conical,
        k_Spherical,
        k_Hemispherical,
        k_Cylindrical,
        k_TaperedCylindrical,
        k_Flame,
        k_InverseConical,
        k_TendFlame
    } shape;
    float baseSize;
    float scale, scaleV, zScale, zScaleV;
    int levels;
    float ratio, ratioPower;
    int lobes;
    float lobeDepth;
    float flare;
    float baseSplits;
    struct stemParameters stemParameters[4];
    struct leafParameters leafParameters;
};

static float
shapeRatio(enum treeParameters::Shape shape, float ratio) {
    switch (shape) {
        case treeParameters::k_Conical:
            return 0.2f + 0.8f * ratio;
            break;
        case treeParameters::k_Spherical:
            return 0.2f + 0.8f * sinf(M_PI * ratio);
            break;
        case treeParameters::k_Hemispherical:
            return 0.2f + 0.8f * sinf(0.5 * M_PI * ratio);
            break;
        case treeParameters::k_Cylindrical:
            return 1.0f;
            break;
        case treeParameters::k_TaperedCylindrical:
            return 0.5f + 0.5f * ratio;
            break;
        case treeParameters::k_Flame:
            return (ratio < 0.7f) ? (ratio / 0.7f) : (1.0f - ratio) / 0.3f;
            break;
        case treeParameters::k_InverseConical:
            return 1.0f - 0.8f * ratio;
            break;
        case treeParameters::k_TendFlame:
            return (ratio < 0.7f) ? (0.5 + 0.5 * ratio / 0.7f) : (0.5 + 0.5f * (1.0 - ratio) / 0.3);
            break;
    }
    return 1;
}


const struct treeParameters CA_BlackOakParameters = {
    treeParameters::k_Hemispherical,
    0.05,
    10, 10, 1, 0,
    3,
    0.018, 1.3,
    5, 0.1,
    1.2,
    2,
    {
        {
            0, 0,
            0, 0, 0,
            1, 0, 0.95,
            0.4, 10, 0,
            8,
            0, 0, 90
        },
        {
            30, -30,
            80, 0, 40,
            0.8, 0.1, 1,
            0.2, 10, 10,
            10,
            40, -70, 150
        },
        {
            45, 10,
            140, 0, 120,
            0.2, 0.05, 1,
            0.1, 10, 10,
            3,
            0, 0, -30
        },
        {
            45, 10,
            140, 0, 5, /* paper has 0, which makes no sense. */
            0.4, 0, 1,
            0, 0, 0,
            1,
            0, 0, 0
        }
    },
    {
        25,
        0.12, 0.66,
        0.8
    }
};

const struct treeParameters BlackTupeloParameters = {
    treeParameters::k_TaperedCylindrical,
    0.2,
    23, 5, 1, 0,
    4,
    0.015, 1.3,
    3, 0.1,
    1,
    0,
    {
        {
            0, 0,
            0, 0, 0,
            1, 0, 1.1,
            0, 0, 0,
            10,
            0, 0, 40
        },
        {
            60, -40,
            140, 0, 50,
            0.3, 0.05, 1,
            0, 0, 0,
            10,
            0, 0, 90
        },
        {
            30, 10,
            140, 0, 25,
            0.6, 0.1, 1,
            0, 0, 0,
            10,
            -10, 0, 150
        },
        {
            45, 10,
            140, 0, 12,
            0.4, 0, 1,
            0, 0, 0,
            1,
            0, 0, 0
        }
    },
    {
        6,
        // The paper has 0.3 here. 30 cm leaves seems way too big - I
        // pulled 0.085 from Wikipedia
        0.085, 0.5,
        0.5
    }
};

const struct treeParameters QuakingAspenParameters = {
    treeParameters::k_TendFlame,
    0.4,
    13, 3, 1, 0,
    3,
    0.015, 1.2,
    5, 0.07,
    0.6,
    0,
    {
        {
            0, 0,
            0, 0, 0,
            1, 0, 1,
            0, 0, 0,
            3,
            0, 0, 20
        },
        {
            60, -50,
            140, 0, 50,
            0.3, 0, 1,
            0, 0, 0,
            5,
            -40, 0, 50
        },
        {
            45, 10,
            140, 0, 30,
            0.6, 0, 1,
            0, 0, 0,
            3,
            -40, 0, 75
        },
        {
            45, 10,
            77, 0, 10,
            /*0*/ 0.4, 0, 1, // paper has 0 for 3Length - surely a typo
            0, 0, 0,
            1,
            0, 0, 0
        }
    },
    {
        25,
//        0.17, 1, // paper has 0.17. Possibly true for young trees, but not true for older aspens
        0.09, 1,
        0.5
    }
};


struct stemData {
    stemData() :
        nf(0) {
    }
        
    int nf;
    std::vector<int> nverts;
    std::vector<int> verts;
    std::vector<float> p;
    std::vector<int> corners;
};

static void
extendstem(const struct arbreParameters &params,
    const struct stemParameters &stem,    
    int level, int segment,
    float baseradius, float tipradius,
    float length, float parentLength, float parentV,
    float branchPower, float previousBranchAngle,
    RtPoint base, RtMatrix coordSystem,
    std::deque<float> &unspreadAngles,
    std::deque<float> &unspreadAxes,
    int previousRing[4],
    unsigned short xsubi[3],
    struct stemData &result) {

    const unsigned short basexsubi[3] = {
        xsubi[0],
        xsubi[1],
        xsubi[2]
    };
    
    int i;
    float v = segment / (float) stem.curveRes;
    float deltaV = 1.0f / stem.curveRes;
    float segmentBaseRadius = baseradius + v * (tipradius - baseradius);    
    float segmentTipRadius = baseradius + (v + deltaV) * (tipradius - baseradius);

    if (params.tree->flare != 0.0f && level == 0) {
        float y = 1.0f - 8.0f * v;
        if (y < 0.0f) y = 0.0f;
        float flare = params.tree->flare * (powf(100.0f, y) - 1.0f) * 0.01f + 1.0f;
        segmentBaseRadius *= flare;

        y = 1.0f - 8.0f * (v + deltaV);
        if (y < 0.0f) y = 0.0f;
        flare = params.tree->flare * (powf(100.0f, y) - 1) * 0.01f + 1.0f;
        segmentTipRadius *= flare;
    }
    
    float segmentLength = length / (float) stem.curveRes;
    float branchAngle = previousBranchAngle;
    float angle;

    RtMatrix nextCoordSystem;

    
    int childlevel = level + 1;
    if (childlevel > 3) childlevel = 3;
    const struct stemParameters &childStem = params.tree->stemParameters[childlevel];

    // Swizzle random seed for branch/leaf growing
    unsigned short branchxsubi[3] = {
        basexsubi[2] + 37,
        basexsubi[1] + 42,
        basexsubi[0] + 79,
    };

    //////////////////////////////////////////////////
    // Branch growing
    //////////////////////////////////////////////////    
    if (level < params.maxLevel && childStem.branches != 0) {

        // Compute the number of branches to generate for this stem
        float maxBranches = 0;
        float childLengthMax = childStem.length;
        float childLength = 0;

        if (level == 0) {
            // The paper is *really* confusing here. It seems to
            // require the number of branches be a function of the
            // length of the child branches, which is a chicken and
            // egg problem for the way the code is structured
            // here. Solution implemented here is to compute the
            // number of branches at v and at v + deltaV and take the
            // max.
            float ratio = (1 - v) / (1 - params.tree->baseSize);
            childLength = childLengthMax * length * shapeRatio(params.tree->shape, ratio);
            float maxBranches1 = childStem.branches * (0.2 + 0.8 * (childLength / length)) /
                childLengthMax / (float) stem.curveRes;
            ratio = (1 - (v + deltaV)) / (1 - params.tree->baseSize);
            childLength = childLengthMax * length * shapeRatio(params.tree->shape, ratio);
            float maxBranches2 = childStem.branches * (0.2 + 0.8 * (childLength / length)) /
                childLengthMax / (float) stem.curveRes;
            maxBranches = maxBranches1;
            if (maxBranches2 > maxBranches) maxBranches = maxBranches2;
        } else {
            maxBranches = childStem.branches * (1.0 - 0.5 * parentV)
                / ((float) stem.curveRes);
        }
        // printf("branches per segment = %f branchPower = %f\n", maxBranches, branchPower);
 

        // Paper: "Any stem that has been cloned or is, itself, a
        // clone reduces its propensity to form clones by half." I
        // assume this actually meant its propensity to form
        // *branches*.
        maxBranches *= branchPower;
        int nBranches = (int) truncf(maxBranches);
        maxBranches -= nBranches;
        if (maxBranches > 0) {
            if (erand48(branchxsubi) < maxBranches) {
                nBranches++;
            }
        }
            
        for (i = 0; i < nBranches; ++i) {
            // branchv is 0 to 1 for the branches in this stem segment
            float branchV = i / (float) nBranches;
            // branchTotalV is 0 to 1 for the entire stem
            float branchTotalV = v + branchV * deltaV;

            if (level == 0) {
                if (branchTotalV < params.tree->baseSize) {
                    continue;
                } else {
                    float ratio = (1 - branchTotalV) / (1 - params.tree->baseSize);
                    childLength = childLengthMax * length * shapeRatio(params.tree->shape, ratio);
                }
            } else {
                childLength = childLengthMax * (length - 0.6 * branchTotalV * length);
            }
            
            if (params.terminalScale != 0 && childLength > 0.01f) {

                
                float childRadius = baseradius * powf(childLength / length, params.tree->ratioPower);

                // "The maximum radius of a stem is explicitly
                // limited to the radius of the parent at the
                // point from which it was spawned."
                float maxChildRadius = segmentBaseRadius + branchV *
                    (segmentTipRadius - segmentBaseRadius);
                if (childRadius > maxChildRadius) childRadius = maxChildRadius;
            
                RiAttributeBegin();

                RiTranslate(base[0] + branchV * segmentLength * coordSystem[2][0],
                    base[1] + branchV * segmentLength * coordSystem[2][1],
                    base[2] + branchV * segmentLength * coordSystem[2][2]);

                branchAngle += childStem.rotate + (2 * erand48(branchxsubi) - 1) *
                    childStem.rotateV;
                RiRotate(branchAngle, 0, 0, 1);
                
                float downAngle;
                if (childStem.downAngleV < 0.0f) {
                    float ratio = (1 - branchTotalV) / (1 - params.tree->baseSize);
                    downAngle = childStem.downAngle + (2 * erand48(branchxsubi) - 1) *
                        childStem.downAngleV * (1 - 2 * shapeRatio(treeParameters::k_Conical, ratio));
                } else {
                    downAngle = childStem.downAngle + (2 * erand48(branchxsubi) - 1) *
                        childStem.downAngleV;
                }
                RiRotate(downAngle, 1, 0, 0);

                if (params.terminalScale != 1.0f &&
                    (level + 1 == params.maxLevel)) {
                    RiScale(params.terminalScale, params.terminalScale, params.terminalScale);
                }
                
                RiConcatTransform(coordSystem);
                
                unsigned short newxsubi[3] = {
                    basexsubi[2] + 37 + i,
                    basexsubi[1] + 42 - i,
                    basexsubi[0] + 79
                };

                if (params.instanceLevel && level + 1 >= params.instanceLevel) {
                    // Find an instanced branch of appropriate size
                    float ratio = childLength / params.maxStemLengths[params.instanceLevel];
                    int bucket = trunc(ratio * 10);
                    int index = bucket * 10 + (erand48(newxsubi) * 10);
                    char buf[256];
                    sprintf(buf, "branch_%d_%d", level + 1, index);
                    RiObjectInstanceV((RtObjectHandle) buf, 0, NULL, NULL);
                } else {
                    growstem(params, childRadius, childLength, level + 1,
                        branchTotalV, newxsubi);
                }
                RiAttributeEnd();
            }
        }
    }

    //////////////////////////////////////////////////
    // Leaf growing
    //////////////////////////////////////////////////
    if (params.leaves && level != 0) {/* && level == params.tree->levels*/ // uncomment to grow only on terminal branches
        const struct leafParameters &leaf = params.tree->leafParameters;
        float maxLeaves = leaf.leaves * shapeRatio(treeParameters::k_TaperedCylindrical, parentV) /
            (float) stem.curveRes;
        int nLeaves = trunc(maxLeaves);
        maxLeaves -= nLeaves;
        if (maxLeaves > 0) {
            if (erand48(branchxsubi) < maxLeaves) {
                nLeaves++;
            }
        }
        for (i = 0; i < nLeaves; ++i) {
            // leafV is 0 to 1 for the leaves in this stem segment
            float leafV = i / (float) nLeaves;

            RiAttributeBegin();
            RiTranslate(base[0] + leafV * segmentLength * coordSystem[2][0],
                base[1] + leafV * segmentLength * coordSystem[2][1],
                base[2] + leafV * segmentLength * coordSystem[2][2]);

            branchAngle += childStem.rotate + (2 * erand48(branchxsubi) - 1) *
                childStem.rotateV;
            RiRotate(branchAngle, 0, 0, 1);            

            if (params.terminalScale != 1.0f &&
                (level + 1 == params.maxLevel)) {
                RiScale(params.terminalScale, params.terminalScale, params.terminalScale);
            }
            
            // I don't think we ever get negative downAngleV for leaves
            float downAngle = childStem.downAngle + (2 * erand48(branchxsubi) - 1) * childStem.downAngleV;
            RiRotate(downAngle, 1, 0, 0);            
            
            RiConcatTransform(coordSystem);

            // Assumes the leaf is modelled with length in Z, width in
            // X
            RiScale(leaf.leafScaleX * leaf.leafScale,
                1.0f,
                leaf.leafScale);

            if (params.instanceLeaves) {
                RiObjectInstanceV((RtObjectHandle) "leaf", 0, NULL, NULL);
            } else {
                RiReadArchive("leaf", NULL, NULL);
            }
            RiAttributeEnd();
        }
    }

    //////////////////////////////////////////////////
    // Stem growing
    //////////////////////////////////////////////////
    
    // Decide whether to split or continue
    int nSplits = trunc(stem.segSplits);
    float splitf = stem.segSplits - nSplits;
    if (splitf > 0) {
        if (erand48(xsubi) < splitf) {
            nSplits++;
        }
    }
    
    // Compute the next segment's coordinate system
    if (segment != stem.curveRes - 1) {
        memcpy(nextCoordSystem, coordSystem, sizeof(RtMatrix));

        // Apply curve variation about local x axis
        if (stem.curveBack == 0) {
            angle = stem.curve / (float) stem.curveRes;
        } else {
            if (segment < stem.curveRes / 2) {
                angle = stem.curve / (0.5f * stem.curveRes);
            } else {
                angle = stem.curveBack / (0.5f * stem.curveRes);
            }
        }
        angle += ((2 * erand48(xsubi) - 1) * stem.curveV) / (float) stem.curveRes;

        rotateMatrix(angle, coordSystem[0][0], coordSystem[0][1],
            coordSystem[0][2], nextCoordSystem);
                

        // Apply spread angle compensations back towards parent(s) stem axes
        for (int j = 0; j < (int) unspreadAngles.size(); ++j) {
            rotateMatrix(unspreadAngles[j], unspreadAxes[j * 3],
                unspreadAxes[j * 3 + 1], unspreadAxes[j * 3 + 2],
                nextCoordSystem);
        }
    }    

    // FIXME: the paper describes ternary splits (nSplits > 1), this
    // only implements binary splits
    if (nSplits && segment != stem.curveRes - 1) {
        // Compute the coordinate systems for the clone and the
        // continuation of the original stem. We rotate the splits
        // around the x axis to minimize tearing.

        // Angle of declination is acosf(local z axis . tree z axis).
        // Since local tree z axis is 0, 0, 1, the dot product is just
        // coordsys[2][2]
        float declination = acosf(coordSystem[2][2]) * 57.2957795131f;

        RtMatrix cloneCoordSystem;
        memcpy(cloneCoordSystem, nextCoordSystem, sizeof(RtMatrix));
        
        std::deque<float> cloneUnspreadAngles(unspreadAngles);
        std::deque<float> cloneUnspreadAxes(unspreadAxes);

        angle = stem.splitAngle - declination;
        if (angle < 0) {
            angle = 0;
        }
        rotateMatrix(angle, coordSystem[0][0], coordSystem[0][1],
            coordSystem[0][2], cloneCoordSystem);
        cloneUnspreadAngles.push_front(-angle / (float) (stem.curveRes - segment + 1));
        cloneUnspreadAxes.push_back(coordSystem[0][0]);
        cloneUnspreadAxes.push_back(coordSystem[0][1]);
        cloneUnspreadAxes.push_back(coordSystem[0][2]);

        rotateMatrix(-angle, coordSystem[0][0], coordSystem[0][1],
            coordSystem[0][2], nextCoordSystem);
        unspreadAngles.push_front(angle / (float) (stem.curveRes - segment + 1));
        unspreadAxes.push_back(coordSystem[0][0]);
        unspreadAxes.push_back(coordSystem[0][1]);
        unspreadAxes.push_back(coordSystem[0][2]);

        // Additional rotation around the tree axis (which is defined
        // to be 0, 0, 1)
        float xi = erand48(xsubi);
        xi *= xi;
        angle = 20 + 0.75 * (30 + fabsf(declination - 90)) * xi;
        if (erand48(xsubi) > 0.5) {
            angle = -angle;
        }
        rotateMatrix(angle, 0, 0, 1, nextCoordSystem);
        int vertbase = result.p.size() / 3;
        
        // Generate split faces
        int newRing[6] = {vertbase, vertbase + 1, vertbase + 2,
                          vertbase + 3, vertbase + 4, vertbase + 5};        

        // Move to end of segment in local Z coordinate system
        RtPoint clonebase;
        clonebase[0] = base[0] + segmentLength * coordSystem[2][0];
        clonebase[1] = base[1] + segmentLength * coordSystem[2][1];
        clonebase[2] = base[2] + segmentLength * coordSystem[2][2];
        base[0] = base[0] + segmentLength * coordSystem[2][0];
        base[1] = base[1] + segmentLength * coordSystem[2][1];
        base[2] = base[2] + segmentLength * coordSystem[2][2];

        for (i = 0; i < 3; ++i) {
            result.p.push_back(clonebase[i] -
                segmentTipRadius * 0.5 * (cloneCoordSystem[0][i] + coordSystem[0][i]) +
                segmentTipRadius * 0.5 * (cloneCoordSystem[1][i] + coordSystem[1][i]));
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(0.5 * (clonebase[i] +
                    segmentTipRadius * 0.5 * (cloneCoordSystem[0][i] + coordSystem[0][i]) +
                    segmentTipRadius * 0.5 * (cloneCoordSystem[1][i] + coordSystem[1][i])) +
                0.5 * (base[i] -
                    segmentTipRadius * 0.5 * (nextCoordSystem[0][i] + coordSystem[0][i]) +
                    segmentTipRadius * 0.5 * (nextCoordSystem[1][i] + coordSystem[1][i])));
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(base[i] +
                segmentTipRadius * 0.5 * (nextCoordSystem[0][i] + coordSystem[0][i]) +
                segmentTipRadius * 0.5 * (nextCoordSystem[1][i] + coordSystem[1][i]));
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(base[i] +
                segmentTipRadius * 0.5 * (nextCoordSystem[0][i] + coordSystem[0][i]) -
                segmentTipRadius * 0.5 * (nextCoordSystem[1][i] + coordSystem[1][i]));            
        }

        for (i = 0; i < 3; ++i) {
            result.p.push_back(0.5 * (clonebase[i] +
                    segmentTipRadius * 0.5 * (cloneCoordSystem[0][i] + coordSystem[0][i]) -
                    segmentTipRadius * 0.5 * (cloneCoordSystem[1][i] + coordSystem[1][i])) +
                0.5 * (base[i] -
                    segmentTipRadius * 0.5 * (nextCoordSystem[0][i] + coordSystem[0][i]) -
                    segmentTipRadius * 0.5 * (nextCoordSystem[1][i] + coordSystem[1][i])));
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(clonebase[i] -
                segmentTipRadius * 0.5 * (cloneCoordSystem[0][i] + coordSystem[0][i]) -
                segmentTipRadius * 0.5 * (cloneCoordSystem[1][i] + coordSystem[1][i]));
        }

        // Connect the dots
        result.nf += 4;
        result.nverts.push_back(5);
        result.nverts.push_back(4);
        result.nverts.push_back(5);
        result.nverts.push_back(4);

        result.verts.push_back(previousRing[0]);
        result.verts.push_back(newRing[0]);
        result.verts.push_back(newRing[1]);
        result.verts.push_back(newRing[2]);
        result.verts.push_back(previousRing[1]);

        result.verts.push_back(previousRing[1]);
        result.verts.push_back(newRing[2]);
        result.verts.push_back(newRing[3]);
        result.verts.push_back(previousRing[2]);

        result.verts.push_back(previousRing[2]);
        result.verts.push_back(newRing[3]);
        result.verts.push_back(newRing[4]);
        result.verts.push_back(newRing[5]);
        result.verts.push_back(previousRing[3]);
        
        result.verts.push_back(previousRing[3]);
        result.verts.push_back(newRing[5]);
        result.verts.push_back(newRing[0]);
        result.verts.push_back(previousRing[0]);


        result.corners.push_back(newRing[1]);
        result.corners.push_back(newRing[4]);
        
        // And continue
        int newRing1[4] = {vertbase, vertbase + 1,
                           vertbase + 4, vertbase + 5};
        int newRing2[4] = {vertbase + 1, vertbase + 2,
                           vertbase + 3, vertbase + 4};
        unsigned short clonexsubi[3] = {
            basexsubi[1] + segment + 137,
            basexsubi[2] - segment + 259,
            basexsubi[0] + 511
        };
        extendstem(params, stem, level, segment + 1, baseradius, tipradius,
            length, parentLength, v + deltaV, branchPower * 0.5f, branchAngle, clonebase, cloneCoordSystem,
            cloneUnspreadAngles, cloneUnspreadAxes, newRing1, clonexsubi, result);
        unsigned short newxsubi[3] = {
            basexsubi[2] + segment + 389,
            basexsubi[1] - segment + 107,
            basexsubi[0] + 989
        };
        extendstem(params, stem, level, segment + 1, baseradius, tipradius,
            length, parentLength, v + deltaV, branchPower * 0.5f, branchAngle, base, nextCoordSystem,
            unspreadAngles, unspreadAxes, newRing2, newxsubi, result);
    } else {
        // Move to end of segment in local Z coordinate system
        base[0] += segmentLength * coordSystem[2][0];
        base[1] += segmentLength * coordSystem[2][1];
        base[2] += segmentLength * coordSystem[2][2];
            
        int vertbase = result.p.size() / 3;

        for (i = 0; i < 3; ++i) {
            result.p.push_back(base[i] - segmentTipRadius * coordSystem[0][i] +
                segmentTipRadius * coordSystem[1][i]);
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(base[i] + segmentTipRadius * coordSystem[0][i] +
                segmentTipRadius * coordSystem[1][i]);
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(base[i] + segmentTipRadius * coordSystem[0][i] -
                segmentTipRadius * coordSystem[1][i]);
        }
        for (i = 0; i < 3; ++i) {
            result.p.push_back(base[i] - segmentTipRadius * coordSystem[0][i] -
                segmentTipRadius * coordSystem[1][i]);
        }

        int newRing[4] = {vertbase, vertbase + 1, vertbase + 2, vertbase + 3};

        // We cap the last segment with one more face
        if (segment == stem.curveRes - 1) {
            result.nf += 5;
        } else {
            result.nf += 4;
        }
        result.nverts.push_back(4);
        result.nverts.push_back(4);
        result.nverts.push_back(4);
        result.nverts.push_back(4);
        if (segment == stem.curveRes - 1) {
            result.nverts.push_back(4);
        }

        result.verts.push_back(previousRing[0]);
        result.verts.push_back(newRing[0]);
        result.verts.push_back(newRing[1]);
        result.verts.push_back(previousRing[1]);

        result.verts.push_back(previousRing[1]);
        result.verts.push_back(newRing[1]);
        result.verts.push_back(newRing[2]);
        result.verts.push_back(previousRing[2]);

        result.verts.push_back(previousRing[2]);
        result.verts.push_back(newRing[2]);
        result.verts.push_back(newRing[3]);
        result.verts.push_back(previousRing[3]);

        result.verts.push_back(previousRing[3]);
        result.verts.push_back(newRing[3]);
        result.verts.push_back(newRing[0]);
        result.verts.push_back(previousRing[0]);

        if (segment == stem.curveRes - 1) {
            result.verts.push_back(newRing[0]);
            result.verts.push_back(newRing[3]);
            result.verts.push_back(newRing[2]);
            result.verts.push_back(newRing[1]);
        } else {
            unsigned short newxsubi[3] = {
                basexsubi[2] + segment + 389,
                basexsubi[1] - segment + 107,
                basexsubi[0] + 989
            };
            extendstem(params, stem, level, segment + 1,
                baseradius, tipradius, length, parentLength, v + deltaV, branchPower,
                branchAngle,
                base, nextCoordSystem,
                unspreadAngles, unspreadAxes, newRing, newxsubi, result);
        }
    }
}

static void
growstem(const struct arbreParameters &params,
    float baseradius, float length, int level, float parentV,
    unsigned short xsubi[3]) {

    const struct stemParameters &stem = params.tree->stemParameters[(level > 3) ? 3 : level];

    int i;
    RtPoint base = { 0, 0, 0 };
    RtMatrix coordSystem = {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 }
    };
    std::deque<float> unspreadAngles;
    std::deque<float> unspreadAxes;
    struct stemData result;
    // Kick things off by adding a ring at the base
    float stemtipradius = params.tree->ratio * params.tree->scale * (1 - stem.taper);
    float radius = baseradius;

    if (level == 0 && params.tree->flare != 0.0f) {
        // Derived from paper equation with z = 0
        radius *= (1 + 0.99 * params.tree->flare); 
    }
    result.p.push_back(-radius);
    result.p.push_back(radius);
    result.p.push_back(0);

    result.p.push_back(radius);
    result.p.push_back(radius);
    result.p.push_back(0);

    result.p.push_back(radius);
    result.p.push_back(-radius);
    result.p.push_back(0);

    result.p.push_back(-radius);
    result.p.push_back(-radius);
    result.p.push_back(0);

    int basering[4] = {0, 1, 2, 3};

    
    RiAttributeBegin();
    extendstem(params, stem, level, 0,
        baseradius, stemtipradius, length, -1 /* no parentLength */,
        parentV, 1.0f, 0.0f, base, coordSystem, unspreadAngles,
        unspreadAxes, basering, xsubi, result);

    RiReadArchive("barkmaterial", NULL, NULL);
    
    RtToken toks[1] = {RI_P};
    RtPointer vals[1] = {&result.p[0]};
    int ncorners = result.corners.size();
    int ncreases = ncorners / 2;
    RtToken *tags = (RtToken *) malloc((ncreases + 2) * sizeof(RtToken));
    RtInt *nargs = (RtInt *) malloc((ncreases + 2) * 2 * sizeof(RtInt));
    RtInt *intargs = (RtInt *) malloc((ncorners + ncreases * 2) * sizeof(RtInt));
    RtFloat *floatargs = (RtFloat *) malloc((1 + ncreases) * sizeof(RtFloat));

    RtToken *tagsptr = tags;
    RtInt *nargsptr = nargs;
    RtInt *intargsptr = intargs;
    RtFloat *floatargsptr = floatargs;
    *tagsptr++ = "interpolateboundary";
    *nargsptr++ = 0;
    *nargsptr++ = 0;

    *tagsptr++ = "corner";
    *nargsptr++ = ncorners;
    *nargsptr++ = 1;
    for (i = 0; i < ncorners; ++i) {
        *intargsptr++ = result.corners[i];
    }
    *floatargsptr++ = 2.0;

    for (i = 0; i < ncreases; ++i) {    
        *tagsptr++ = "crease";
        *nargsptr++ = 2;
        *nargsptr++ = 1;
        *intargsptr++ = result.corners[2 * i];
        *intargsptr++ = result.corners[2 * i + 1];
        *floatargsptr++ = 5.0;
    }

    if (result.nf) {
#if 0
        // Useful for debugging subdivs
        RiPointsPolygonsV(result.nf, &result.nverts[0],
            &result.verts[0],
            1, toks, vals);
#else
        RiSubdivisionMeshV(RI_CATMARK, result.nf, &result.nverts[0],
            &result.verts[0],
            2 + ncreases, tags, nargs, intargs, floatargs,
            1, toks, vals);
#endif
    }
    free(tags);
    free(nargs);
    free(intargs);
    free(floatargs);
    RiAttributeEnd();
}

PRMANEXPORT RtVoid Subdivide2(RtContextHandle ctx, RtFloat detail,
    RtInt n, RtToken toks[], RtPointer vals[]) {

    float length = -1, radius, v = 0.0f;
    int i, j;
    struct arbreParameters params;
    RtToken type = "tupelo";
    RtInt seed = 42; // Arbitrary default seed

    params.maxLevel = -1;
    params.startLevel = 0;
    params.instanceLevel = -1;
    params.leaves = 1;
    params.instanceLeaves = 1;
    params.terminalScale = 1.0;
    params.tree = &BlackTupeloParameters;
    for (i = 0; i < n; i++)
    {
        if (strstr(toks[i], "Seed")) {
            seed = *(RtInt *) vals[i];
        }
        else if (strstr(toks[i], "MaxLevel")) {
            params.maxLevel = *(RtInt *) vals[i];
        }
        else if (strstr(toks[i], "StartLevel")) {
            params.startLevel = *(RtInt *) vals[i];
        }
        else if (strstr(toks[i], "InstanceLevel")) {
            params.instanceLevel = *(RtInt *) vals[i];
        }
        else if (strstr(toks[i], "InstanceLeaves")) {
            params.instanceLeaves = *(RtInt *) vals[i];
        }
        else if (strstr(toks[i], "Leaves")) {
            params.leaves = *(RtInt *) vals[i];
        }
        else if (strstr(toks[i], "TerminalScale")) {
            params.terminalScale = *(RtFloat *) vals[i];
        }
        else if (strstr(toks[i], "Type")) {
            type = *(RtToken *) vals[i];
            if (strstr(type, "oak")) {
                params.tree = &CA_BlackOakParameters;
            } else if (strstr(type, "tupelo")) {
                params.tree = &BlackTupeloParameters;
            } else if (strstr(type, "aspen")) {
                params.tree = &QuakingAspenParameters;
            }
            
        }
        else if (strstr(toks[i], "LengthOverride")) {
            length = *(RtFloat *) vals[i];
        }
    }

    if (params.maxLevel == -1) {
        params.maxLevel = params.tree->levels;
    } else if (params.maxLevel > params.tree->levels) {
        params.maxLevel = params.tree->levels;
    }

    if (params.instanceLevel == -1) {
        // By default, instance the last level
        params.instanceLevel = params.tree->levels;
    } else if (params.instanceLevel > params.tree->levels) {
        params.instanceLevel = params.tree->levels;
    }
    
    if (params.startLevel < 0) {
        params.startLevel = 0;
    } else if (params.startLevel >= params.tree->levels) {
        params.startLevel = params.tree->levels;
    }

    // Compute maximum branch lengths at each level
    params.maxStemLengths[0] = params.tree->scale + params.tree->scaleV;
    for (i = 1; i <= params.tree->levels; ++i) {
        params.maxStemLengths[i] = params.maxStemLengths[i - 1] *
            params.tree->stemParameters[i > 3 ? 3 : i].length;
    }
    unsigned short xsubi[3];

    if (params.instanceLevel != 0) {
        char buf[256];
        RixContext *rix = RixGetContext();
        RixStorage* storage = (RixStorage*) rix->GetRixInterface(k_RixGlobalData);
        storage->Lock();
        sprintf(buf, "arbre%s", type);
        void* init = storage->Get(buf);
        if (init == NULL) {
            printf("  doing once-only initialization of shared branches..\n");
        
            for (int instanceLevel = params.tree->levels;
                  instanceLevel >= params.instanceLevel;
                 --instanceLevel) {

                printf("    initing at level %d", instanceLevel);

                // The radius of a branch grown at a particular level is a
                // function of where it is grown along the parent, but is
                // also dependent on where its parent is grown relative to
                // its grandparent. We won't be able to correctly deal
                // with the parent offset relative to the grandparent with
                // instancing. We will approximate it by using an average:
                // pretending that the parent was grown halfway along the
                // grandparent.
                float childLength, childLengthMax, childRadius;
                float parentLength = params.tree->scale;
                float parentRadius = params.maxStemLengths[0] * params.tree->ratio;
                float parentLengthMax;
                for (j = 1; j < instanceLevel; ++j) {
                    childLengthMax = params.tree->stemParameters[(j > 3) ? 3 : j].length;

                    if (j == 1) {
                        float ratio = (1 - 0.5f) / (1 - params.tree->baseSize);                    
                        childLength = childLengthMax * parentLength * shapeRatio(params.tree->shape, ratio);
                    } else {
                        childLength = childLengthMax * (parentLength - 0.6f * 0.5f * parentLength);
                    }
                    childRadius = parentRadius * powf(childLength / parentLength,
                        params.tree->ratioPower);
                    parentRadius = childRadius;
                    parentLength = childLength;
                }

                parentLengthMax = params.maxStemLengths[instanceLevel - 1];
                childLengthMax = params.maxStemLengths[instanceLevel];
            
                // Create instanced branches of various lengths
                for (i = 0; i < 100; i++) {
                    xsubi[0] = i;
                    xsubi[1] = i;
                    xsubi[2] = i;
        
                    float bucket = trunc((i + 10.0f) / 10.0f);
                    RtToken newtoks[1] = {"string __handleid"};
                    RtPointer newvals[1];
                    RtToken nameptr = buf;
                    sprintf(buf, "branch_%d_%d", instanceLevel, i);
                    newvals[0] = &nameptr;
                    RiObjectBeginV(1, newtoks, newvals);

                    childLength = bucket * 0.1f * childLengthMax;
                    childRadius = parentRadius * powf(childLength / parentLength, params.tree->ratioPower);

                    // FIXME: this is a hack. We need to pass a parentV to
                    // growstem, because it controls the branching and
                    // leaf density of grandchildren. We could compute it
                    // "correctly" here, but it would often be out of
                    // range given the shenanigans already going on with
                    // the max child/parent. Computing it like this is a
                    // compromise that seems to produce acceptable
                    // results.
                    float parentV = 1.0f - bucket * 0.1f;
                    growstem(params, childRadius, childLength, instanceLevel, parentV, xsubi);
                    RiObjectEnd();
                }
                init = malloc(sizeof(int));
                sprintf(buf, "arbre%s", type);                
                storage->Set(buf, init, NULL);
                printf(".. done!\n");
                fflush(stdout);
            }
        }
        storage->Unlock();
    }
    
    xsubi[0] = seed;
    xsubi[1] = seed;
    xsubi[2] = seed;

    // Compute length if it hasn't already been set
    if (length < 0) {
        if (params.startLevel == 0) {
            length = params.tree->scale + (2 * erand48(xsubi) - 1) * params.tree->scaleV;
            radius = params.tree->ratio * length;
        } else {
            // Similar hackery as before
            float childLength, childRadius;
            float parentLength = params.tree->scale;
            float parentRadius = params.maxStemLengths[0] * params.tree->ratio;
            float childLengthMax = params.tree->stemParameters[1].length;
            for (j = 1; j < params.startLevel; ++j) {
                childLengthMax = params.tree->stemParameters[(j > 3) ? 3 : j].length;

                if (j == 1) {
                    float ratio = (1 - 0.5f) / (1 - params.tree->baseSize);                    
                    childLength = childLengthMax * parentLength * shapeRatio(params.tree->shape, ratio);
                } else {
                    childLength = childLengthMax * (parentLength - 0.6f * 0.5f * parentLength);
                }
                childRadius = parentRadius * powf(childLength / parentLength,
                    params.tree->ratioPower);
                parentRadius = childRadius;
                parentLength = childLength;
            }
            float xi = erand48(xsubi);
            length = xi * childLengthMax * (parentLength - 0.6f * 0.5f * parentLength);
            radius = parentRadius * powf(length / parentLength,
                    params.tree->ratioPower);
            v = 1 - xi;
        }
    }
    // FIXME: if length is supplied by the user, the v here starts at
    // zero (and shouldn't)
    growstem(params, params.tree->ratio * length, length, params.startLevel, v, xsubi);
}
