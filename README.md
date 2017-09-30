# arbre
This is a partial implementation of Jason Weber and Joseph Penn's paper,
"Creation and Rendering of Realistic Trees", from SIGGRAPH 1995. This code
should not be treated as an example on how to implement that paper properly
(even though there are some confusing sections in the paper that I've tried to
interpret in a useful fashion), but more as a basic example on how one could
set up nested object instance hierarchies to maximize memory efficiency in
PRMan 17.0 (as described in my User's Group Meeting talk at SIGGRAPH 2012).

Sections 4.1 through 4.5 of the Penn/Weber paper are mostly implemented, 4.6
through 4.9 are not. I also did not implement ternary stem splits, handling of
-nCurveV (helical stems), or trunk lobes.

To use this plugin, compile it as a DSO and load it into
RiProcedural2:

Procedural2 "DynamicLoad" "SimpleBound" "float[6] bound" [-20 20
    -20 20 0 20] "string dsoname" ["arbre"] "int Seed" [66] "int
    MaxLevel" [5] "int StartLevel" [0] "int InstanceLevel" [4] "int
    Leaves" [1] "int InstanceLeaves" [1]

The useful parameters are:

"string Type": tree type. Currently supported: oak, tupelo, and aspen.
   The parameters are straight out of the paper, slightly modified
   in a few cases

"int Seed": random seed for the tree. If you are growing a forest
   of trees, every tree should have its own seed

"int MaxLevel": maximum level of the tree - this overrides the Levels
   parameter in the paper

"int InstanceLevel": level at which to use instanced branches. For
   example, if InstanceLevel == MaxLevel, the terminal level of the
   tree will use instanced branches. For two levels of terminal instanced
   branches, set InstanceLevel == MaxLevel - 1, and so on. To disable
   instanced branches, set this parameter to 0.

"int Leaves": toggle switch for whether leaves will be grown or not

"int InstanceLeaves": whether leaves will be instanced. If 1,
   the procedural requires a "leaf" ObjectInstance to have
   been previously defined. If 0, the procedural requires a "leaf"
   inline archive to be previously defined.

"int StartLevel": by default, the procedural starts at Level 0, this
   parameter overrides that default. This is mainly useful for
   visualization purposes


There are several known bugs with this plugin that ought to be addressed
before even considering using this for production work.

- If a stem has high curvature, its children may be detached at the base. This
is "fixed" by switching to polygons. The proper fix for subdivision meshes is
to interpolate using the cubic b-spline rather than just perform linear
interpolation as is done here

- The actual formulae for the leaf/branch frequency/stem length/stem radius
are dependent not only parent properties but also grandparent properties. So,
the simple-minded instancing being used here where we just look for the branch
of the right length isn't close to correct as far the Penn/Weber model is
concerned. I've attempted to compensate for this by associating an arbitrary
parent offset with each stem length, and this seems to be "good enough" at a
first glance. A more rigorous derivation, coupled perhaps by searching for
appropriate instances in multiple dimensions (e.g. radius and length), may
provide more accurate results.

- The use of erand48 throughout is completely ad hoc and wacky.

Naturally, there are other useful features that should be implemented as well
(like switching amongst multiple leaf types).


Julian Fong <jfong@pixar.com>
