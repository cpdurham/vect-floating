module Data.Vect.Floating.Util.Dim3 where

import Data.Vect.Floating.Base

--------------------------------------------------------------------------------

-- | Example: @structVec3 [1,2,3,4,5,6] = [ Vec3 1 2 3 , Vec3 4 5 6]@.
structVec3 :: [a] -> [Vec3 a]
structVec3 [] = []
structVec3 (x:y:z:ls) = (Vec3 x y z):(structVec3 ls) 
structVec3 _ = error "structVec3"

-- | The opposite of "structVec3".
destructVec3 :: [Vec3 a] -> [a]
destructVec3 [] = []
destructVec3 ((Vec3 x y z):ls) = x:y:z:(destructVec3 ls)  

--------------------------------------------------------------------------------

det3 :: Floating a => Vec3 a -> Vec3 a -> Vec3 a -> a
det3 u v w = det (u,v,w)

--------------------------------------------------------------------------------

translate3X :: Num a => a -> Vec3 a -> Vec3 a
translate3Y :: Num a => a -> Vec3 a -> Vec3 a
translate3Z :: Num a => a -> Vec3 a -> Vec3 a

translate3X t (Vec3 x y z) = Vec3 (x+t) y z 
translate3Y t (Vec3 x y z) = Vec3 x (y+t) z 
translate3Z t (Vec3 x y z) = Vec3 x y (z+t) 

vec3X :: Num a => Vec3 a
vec3Y :: Num a => Vec3 a
vec3Z :: Num a => Vec3 a

vec3X = Vec3 1 0 0
vec3Y = Vec3 0 1 0
vec3Z = Vec3 0 0 1

rotMatrixZ :: Floating a => a -> Mat3 a
rotMatrixY :: Floating a => a -> Mat3 a
rotMatrixX :: Floating a => a -> Mat3 a

-- These are intended for multiplication on the /right/.
-- Should be consistent with the rotation around an arbitrary axis 
-- (eg, @rotMatrixY a == rotate3 a vec3Y@)
rotMatrixZ a = Mat3 (Vec3 c s 0) (Vec3 (-s) c 0) (Vec3 0 0 1) where c = cos a; s = sin a
rotMatrixY a = Mat3 (Vec3 c 0 (-s)) (Vec3 0 1 0) (Vec3 s 0 c) where c = cos a; s = sin a
rotMatrixX a = Mat3 (Vec3 1 0 0) (Vec3 0 c s) (Vec3 0 (-s) c) where c = cos a; s = sin a

--------------------------------------------------------------------------------

rotate3' {- ' CPP is sensitive to primes -} 
  :: Floating a => a       -- ^ angle (in radians)
  -> Normal3 a   -- ^ axis (should be a /unit/ vector!) 
  -> Vec3 a      -- ^ vector
  -> Vec3 a      -- ^ result
rotate3' angle axis v = v .* (rotMatrix3' axis angle)

rotate3 
  :: Floating a => a    -- ^ angle (in radians)
  -> Vec3 a   -- ^ axis (arbitrary nonzero vector)
  -> Vec3 a   -- ^ vector
  -> Vec3 a  -- ^ result
rotate3 angle axis v = v .* (rotMatrix3 axis angle)
      
-- | Rotation around an arbitrary 3D vector. The resulting 3x3 matrix is intended for multiplication on the /right/. 
rotMatrix3 :: Floating a => Vec3 a -> a -> Mat3 a
rotMatrix3 v a = rotMatrix3' (mkNormal v) a

rotMatrixOrtho3 :: Floating a => Vec3 a -> a -> Ortho3 a
rotMatrixOrtho3 v a = toOrthoUnsafe $ rotMatrix3 v a

-- | Rotation around an arbitrary 3D /unit/ vector. The resulting 3x3 matrix is intended for multiplication on the /right/. 
rotMatrix3' :: {- ' CPP is sensitive to primes -} Floating a => Normal3 a -> a -> Mat3 a
rotMatrix3' u a = 
  let v = fromNormal u
      c = cos a
      s = sin a
      m1 = scalarMul (1-c) (outer v v)
      x = _1 v
      y = _2 v
      z = _3 v
      m2 = Mat3 (Vec3   c    ( s*z) (-s*y))
                (Vec3 (-s*z)   c    ( s*x))
                (Vec3 ( s*y) (-s*x)   c   )
  in (m1 &+ m2)

rotMatrixOrtho3' :: {- ' CPP is sensitive to primes -} Floating a => Normal3 a -> a -> Ortho3 a
rotMatrixOrtho3' u a = toOrthoUnsafe $ rotMatrix3' u a

--------------------------------------------------------------------------------

-- | Reflects a vector to an axis: that is, the result of @reflect n v@ is
-- 2\<n,v\>n - v
reflect :: Floating a => Normal3 a -> Vec3 a -> Vec3 a
reflect u v = (s *& n) &- v where 
  n = fromNormal u
  s = 2 * (n &. v)

reflect' :: Floating a => Normal3 a -> Normal3 a -> Normal3 a
reflect' u x = toNormalUnsafe $ reflect u (fromNormal x)
  
refract :: (Floating a, Ord a) => a -> Normal3 a -> Vec3 a -> Vec3 a
refract eta u v = s *& fromNormal w where
  s = norm v 
  w = refract' eta u (toNormalUnsafe $ v &* (1.0/s))
  
-- | Refraction.
-- First parameter (@eta@) is the relative refraction index 
--
-- >        refl_inside
-- > eta = --------------
-- >        refl_outside
--
-- where \"inside\" is the direction of the second argument 
-- (to vector normal to plane which models the boundary 
-- between the two materials). That is, total internal reflection
-- can occur when @eta>1@.
--
-- The convention is that the origin is the point of intersection
-- of the ray and the surface, and all the vectors \"point away\"
-- from here (unlike, say, GLSL's @refract@, where the incident
-- vector \"points towards\" the material)
refract' {- ' CPP is sensitive to primes -} 
  :: (Floating a, Ord a) => a -> Normal3 a -> Normal3 a -> Normal3 a
refract' eta u i = 
  if k<0
    then reflect' u i 
    else toNormalUnsafe $ ((-eta) *& v) &- (- eta*c + sqrt k) *& n 
  where
    n = fromNormal u
    v = fromNormal i
    c = n &. v
    k = 1 - eta*eta*(1-c*c)

-- | When total internal reflection would occur, we return "Nothing".
refractOnly' {- ' CPP is sensitive to primes -}  
  :: (Floating a, Ord a) => a -> Normal3 a -> Normal3 a -> Maybe (Normal3 a)
refractOnly' eta u i = 
  if k<0
    then Nothing 
    else Just $ toNormalUnsafe $ ((-eta) *& v) &- (- eta*c + sqrt k) *& n 
  where
    n = fromNormal u
    v = fromNormal i
    c = n &. v
    k = 1 - eta*eta*(1-c*c)

--------------------------------------------------------------------------------
