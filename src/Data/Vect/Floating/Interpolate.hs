{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses #-}

-- TODO: interpolation for Ortho3 matrices using the (short) quaternion 'slerpU'

-- | Interpolation of vectors. 
-- Note: we interpolate unit vectors differently from ordinary vectors.

module Data.Vect.Floating.Interpolate where

--------------------------------------------------------------------------------

import Data.Vect.Floating.Base

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

  
