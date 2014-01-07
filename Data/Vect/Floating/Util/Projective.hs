{-# LANGUAGE FlexibleContexts, ScopedTypeVariables #-}
-- | Classic 4x4 projective matrices, encoding the affine transformations of R^3.
-- Our convention is that they are intended for multiplication on
-- the /right/, that is, they are of the form
--
-- >     _____
-- > [  |     |  0  ]
-- > [  | 3x3 |  0  ]
-- > [  |_____|  0  ]
-- > [  p  q  r  1  ]
--
-- Please note that by default, OpenGL stores the matrices (in memory) by columns, while we 
-- store them by rows; but OpenGL also use the opposite convention (so the OpenGL projective matrices 
-- are intended for multiplication on the /left/). So in effect, they are the same when stored in the memory,
-- say with @poke :: Ptr Mat4 -> Mat4 -> IO ()@.
--
-- Warning: The naming conventions will probably change in the future.

module Data.Vect.Floating.Util.Projective where

import Data.Vect.Floating.Base
import Data.Vect.Floating.Util.Dim3

import qualified Data.Vect.Floating.Util.Dim4 as Dim4

--------------------------------------------------------------------------------

rotMatrixProj4' :: {- ' CPP is sensitive to primes -} (Floating a, Projective a Vec3 Mat3 Ortho3 b Proj4) => a -> Normal3 a -> Proj4 a
rotMatrixProj4' angle axis = linear $ rotMatrix3' axis angle

rotMatrixProj4 :: Floating a => a -> Vec3 a -> Proj4 a
rotMatrixProj4 angle axis = linear $ rotMatrix3 axis angle

-- | synonym for "rotateAfterProj4"
rotateProj4 :: Floating a => a -> Normal3 a -> Proj4 a -> Proj4 a
rotateProj4 = rotateAfterProj4

-- | Synonym for @\m -> m .*. rotMatrixProj4 angle axis@.
rotateAfterProj4 :: Floating a => a -> Normal3 a -> Proj4 a -> Proj4 a
rotateAfterProj4 angle axis m = m .*. (rotMatrixProj4' angle axis) 

-- | Synonym for @\m -> rotMatrixProj4 angle axis .*. m@.
rotateBeforeProj4 :: Floating a => a -> Normal3 a -> Proj4 a -> Proj4 a
rotateBeforeProj4 angle axis m = (rotMatrixProj4' angle axis) .*. m 

---------------

--scalingUniformProj3 :: Flt -> Proj3
--scalingUniformProj3 x = scaling (Vec2 x x)

scalingUniformProj4 :: Floating a => a -> Proj4 a
scalingUniformProj4 x = scaling (Vec3 x x x)

-- | Equivalent to @\m -> scaling v .*. m@.
scaleBeforeProj4 :: Floating a => Vec3 a -> Proj4 a -> Proj4 a
scaleBeforeProj4 (Vec3 u v w) p4 =  
  toProjectiveUnsafe $ 
    Mat4 (u*&a) (v*&b) (w*&c) t
  where
    Mat4 a b c t = fromProjective p4

-- | Equivalent to @\m -> m .*. scaling v@.
scaleAfterProj4 :: Floating a => Vec3 a -> Proj4 a -> Proj4 a
scaleAfterProj4 v p4 =
  toProjectiveUnsafe $ 
    Mat4 (a&!w) (b&!w) (c&!w) (t&!w)
  where
    w = extendWith 1 v
    Mat4 a b c t = fromProjective p4
    
---------------

-- | Synonym for "translateAfter4"
translate4 :: Floating a => Vec3 a -> Proj4 a -> Proj4 a
translate4 = translateAfter4

-- | Equivalent to @\m -> m .*. translation v@.
translateAfter4 :: Floating a => Vec3 a -> Proj4 a -> Proj4 a
translateAfter4 v p4 = 
  toProjectiveUnsafe $
    Mat4 r1 r2 r3 (extendWith 0 v &+ r4)
  where
    Mat4 r1 r2 r3 r4 = fromProjective p4 

-- | Equivalent to @\m -> translation v .*. m@.
translateBefore4 :: (Floating a, Extend a Mat3 Mat4, RightModule (Vec3 a) (Mat3 a)) => Vec3 a -> Proj4 a -> Proj4 a
translateBefore4 v p4 = 
  toProjectiveUnsafe $ 
    Mat4 r1 r2 r3 (extendWith 0 u &+ r4) 
  where 
   u = v .*  trim mat 
   mat@(Mat4 r1 r2 r3 r4) = fromProjective p4
   
---------------

