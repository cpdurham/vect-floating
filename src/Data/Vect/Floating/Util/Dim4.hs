-- | Rotation around an arbitrary plane in four dimensions, and other miscellanea.
-- Not very useful for most people, and not re-exported by "Data.Vect".

module Data.Vect.Floating.Util.Dim4 where

import Data.Vect.Floating.Base
import Data.Vect.Floating.GramSchmidt

structVec4 :: [a] -> [Vec4 a]
structVec4 [] = []
structVec4 (x:y:z:w:ls) = (Vec4 x y z w):(structVec4 ls) 
structVec4 _ = error "structVec4"

destructVec4 :: [Vec4 a] -> [a]
destructVec4 [] = []
destructVec4 ((Vec4 x y z w):ls) = x:y:z:w:(destructVec4 ls)  

--det4 :: Vec4 -> Vec4 -> Vec4 -> Vec4 -> Flt
--det4 u v w z = det (u,v,w,z)

translate4X :: Num a => a -> Vec4 a -> Vec4 a
translate4Y :: Num a => a -> Vec4 a -> Vec4 a
translate4Z :: Num a => a -> Vec4 a -> Vec4 a
translate4W :: Num a => a -> Vec4 a -> Vec4 a

translate4X t (Vec4 x y z w) = Vec4 (x+t) y z w 
translate4Y t (Vec4 x y z w) = Vec4 x (y+t) z w 
translate4Z t (Vec4 x y z w) = Vec4 x y (z+t) w
translate4W t (Vec4 x y z w) = Vec4 x y z (w+t) 

vec4X :: Num a => Vec4 a
vec4Y :: Num a => Vec4 a
vec4Z :: Num a => Vec4 a
vec4W :: Num a => Vec4 a

vec4X = Vec4 1 0 0 0
vec4Y = Vec4 0 1 0 0
vec4Z = Vec4 0 0 1 0
vec4W = Vec4 0 0 0 1

--------------------------------------------------------------------------------

-- |If @(x,y,u,v)@ is an orthonormal system, then (written in pseudo-code)
-- @biVector4 (x,y) = plusMinus (reverse $ biVector4 (u,v))@.
-- This is a helper function for the 4 dimensional rotation code.
-- If @(x,y,z,p,q,r) = biVector4 a b@, then the corresponding antisymmetric tensor is
--
-- > [  0  r  q  p ]
-- > [ -r  0  z -y ]
-- > [ -q -z  0  x ]
-- > [ -p  y -x  0 ]
biVector4 :: Num a => Vec4 a -> Vec4 a -> (a,a,a,a,a,a)
biVector4 (Vec4 x y z w) (Vec4 a b c d) = 
  ( x*b-y*a , x*c-z*a , x*d-w*a , y*c-z*b , -y*d+w*b , z*d-w*c )

-- | the corresponding antisymmetric tensor
biVector4AsTensor :: Num a => Vec4 a -> Vec4 a -> Mat4 a
biVector4AsTensor v w = 
  Mat4 ( Vec4   0  ( r) ( q) ( p) )
       ( Vec4 (-r)   0  ( z) (-y) )
       ( Vec4 (-q) (-z)   0  ( x) )
       ( Vec4 (-p) ( y) (-x)   0  )
  where 
    (x,y,z,p,q,r) = biVector4 v w

-- | We assume that the axes are normalized and /orthogonal/ to each other!
rotate4' :: {- ' CPP is sensitive to primes -} Floating a => a -> (Normal4 a,Normal4 a) -> Vec4 a -> Vec4 a
rotate4' angle axes v = v .* (rotMatrix4' angle axes)

-- | We assume only that the axes are independent vectors.
rotate4 :: Floating a => a -> (Vec4 a,Vec4 a) -> Vec4 a -> Vec4 a
rotate4 angle axes v = v .* (rotMatrix4 angle axes)

-- | Rotation matrix around a plane specified by two normalized and /orthogonal/ vectors.
-- Intended for multiplication on the /right/!
rotMatrix4' :: {- ' CPP is sensitive to primes -} Floating a => a -> (Normal4 a,Normal4 a) -> Mat4 a
rotMatrix4' angle (u1,u2) = m1 &+ (s *& m2) &+ m3 
  where
    v = fromNormal u1 ; w = fromNormal u2
    c = cos angle ; s = sin angle
    m1 = scalarMul (1-c) ( outer v v  &+  outer w w )
    m2 = biVector4AsTensor v w
    m3 = diag (Vec4 c c c c)

-- | We assume only that the axes are independent vectors.
rotMatrix4 :: Floating a => a -> (Vec4 a,Vec4 a) -> Mat4 a  
rotMatrix4 angle axes = 
  rotMatrix4' angle $ liftPair toNormalUnsafe $ gramSchmidtNormalize axes 
  where 
    liftPair f (x,y) = (f x, f y)
    
    
