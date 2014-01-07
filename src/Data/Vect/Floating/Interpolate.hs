{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}

-- TODO: interpolation for Ortho3 matrices using the (short) quaternion 'slerpU'

-- | Interpolation of vectors. 
-- Note: we interpolate unit vectors differently from ordinary vectors.

module Data.Vect.Floating.Interpolate where

--------------------------------------------------------------------------------

import Data.Vect.Floating.Base
import Data.Vect.Floating.Util.Dim2 (sinCos',angle2')
import Data.Vect.Floating.Util.Dim3 (rotate3')

--------------------------------------------------------------------------------

class Interpolate a v where
  interpolate :: a -> v -> v -> v
  
instance Num a => Interpolate a a where
  interpolate t x y = x + t*(y-x)

--------------------------------------------------------------------------------

instance Floating a => Interpolate a (Vec2 a) where interpolate t x y = x &+ t *& (y &- x)
instance Floating a => Interpolate a (Vec3 a) where interpolate t x y = x &+ t *& (y &- x)
instance Floating a => Interpolate a (Vec4 a) where interpolate t x y = x &+ t *& (y &- x)

--------------------------------------------------------------------------------

{-
instance Interpolate Normal2 where
  interpolate t nx ny = sinCos' $ ax + t*adiff where
    ax = angle2' nx
    ay = angle2' ny
    adiff = helper (ay - ax)
    helper d 
      | d < -pi   = d + twopi
      | d >  pi   = d - twopi
      | otherwise = d
    twopi = 2*pi
    
instance Interpolate Normal3 where 
  interpolate t nx ny = 
    if maxAngle < 0.001  -- more or less ad-hoc critical angle
      then mkNormal $ interpolate t x y
      else toNormalUnsafe $ rotate3' (t*maxAngle) (mkNormal axis) x where
    x = fromNormal nx
    y = fromNormal ny
    axis = (x &^ y)
    maxAngle = acos (x &. y)
-}        

instance Floating a => Interpolate a (Normal2 a) where interpolate = slerp
instance Floating a => Interpolate a (Normal3 a) where interpolate = slerp
instance Floating a => Interpolate a (Normal4 a) where interpolate = slerp
        
--------------------------------------------------------------------------------
  
{-# SPECIALIZE slerp :: Float -> Normal2 Float -> Normal2 Float -> Normal2 Float #-}
{-# SPECIALIZE slerp :: Float -> Normal3 Float -> Normal3 Float -> Normal3 Float #-}
{-# SPECIALIZE slerp :: Float -> Normal4 Float -> Normal4 Float -> Normal4 Float #-}
    
{-# SPECIALIZE slerp :: Double -> Normal2 Double -> Normal2 Double -> Normal2 Double #-}
{-# SPECIALIZE slerp :: Double -> Normal3 Double -> Normal3 Double -> Normal3 Double #-}
{-# SPECIALIZE slerp :: Double -> Normal4 Double -> Normal4 Double -> Normal4 Double #-}

-- | Spherical linear interpolation.
-- See <http://en.wikipedia.org/wiki/Slerp>    
slerp :: UnitVector a v u => a -> u a -> u a -> u a
slerp t n0 n1 = toNormalUnsafe v where
  v = (p0 &* y0) &+ (p1 &* y1) 
  p0 = fromNormal n0
  p1 = fromNormal n1
  omega = acos (p0 &. p1)
  s = sin omega
  y0 = sin (omega*(1-t)) / s 
  y1 = sin (omega*   t ) / s
  
--------------------------------------------------------------------------------

  
