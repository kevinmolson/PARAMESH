#if 0
  The first three constants in this group are used to specify the edges of 
  cells or blocks, where for cell they have the obvious meaning. When used
  in connection with blocks, LEFT_EDGE 
  indicates guard cells to the left of the block, RIGHT_EDGE indicates
  guard cells to the right of the block and CENTER indicates the interior 
  of the block.
#endif

#define LEFT_EDGE 1
#define CENTER 2
#define RIGHT_EDGE 3
#define WHOLE_VECTOR 4
#define NO_VEC 5
#define ALLVARS -1


#if 0
  MDIM defines the maximum number of dimensions
  supported in the code.
#endif

#define MDIM 3


#if 0
  The next two constants are used in connection of integer block boundaries
  LOW refers to the lowest index and HIGH refers to the highest index of the 
  block. These can be used interchangeably for either the whole block 
  including guardcells or for the interior only.
#endif

#define LOW 1
#define HIGH 2

