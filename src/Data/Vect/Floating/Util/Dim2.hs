module Data.Vect.Floating.Util.Dim2 where

import Data.Vect.Floating.Base

-- | Example: @structVec2 [1,2,3,4] = [ Vec2 1 2 , Vec2 3 4 ]@.
structVec2 :: [a] -> [Vec2 a]
structVec2 [] = []
structVec2 (x:y:ls) = (Vec2 x y):(structVec2 ls) 
structVec2 _ = error "structVec2"

-- | The opposite of "structVec2".
destructVec2 :: [Vec2 a] -> [a]
destructVec2 [] = []
destructVec2 ((Vec2 x y):ls) = x:y:(destructVec2 ls)  

det2 :: Floating a => Vec2 a -> Vec2 a -> a
det2 u v = det (u,v)

vec2X :: Num a => Vec2 a
vec2Y :: Num a => Vec2 a

vec2X = Vec2 1 0 
vec2Y = Vec2 0 1 

translate2X :: Num a => a -> Vec2 a -> Vec2 a
translate2Y :: Num a => a -> Vec2 a -> Vec2 a

translate2X t (Vec2 x y) = Vec2 (x+t) y 
translate2Y t (Vec2 x y) = Vec2 x (y+t) 

-- | unit vector with given angle relative to the positive X axis (in the positive direction, that is, CCW).
-- A more precise name would be @cosSin@, but that sounds bad :)
sinCos :: Floating a => a -> Vec2 a
sinCos a = Vec2 (cos a) (sin a)

sinCos' {- ' CPP is sensitive to primes -} :: Floating a => a -> Normal2 a
sinCos' = toNormalUnsafe . sinCos

sinCosRadius :: Floating a => a    -- ^ angle (in radians)
             -> a    -- ^ radius
             -> Vec2 a
sinCosRadius a r = Vec2 (r * cos a) (r * sin a)

-- | The angle relative to the positive X axis
angle2 :: RealFloat a => Vec2 a -> a
angle2 (Vec2 x y) = atan2 y x

angle2' {- ' CPP is sensitive to primes -} :: RealFloat a => Normal2 a -> a
angle2' = angle2 . fromNormal

-- | Rotation matrix by a given angle (in radians), counterclockwise.
rotMatrix2 :: Floating a => a -> Mat2 a
rotMatrix2 a = Mat2 (Vec2 c s) (Vec2 (-s) c) where c = cos a; s = sin a

rotMatrixOrtho2 :: Floating a => a -> Ortho2 a
rotMatrixOrtho2 = toOrthoUnsafe . rotMatrix2

rotate2 :: Floating a => a -> Vec2 a -> Vec2 a
rotate2 a v = v .* (rotMatrix2 a) 

-- |Rotates counterclockwise by 90 degrees.
rotateCCW :: Floating a => Vec2 a -> Vec2 a
rotateCCW (Vec2 x y) = Vec2 (-y) x

-- |Rotates clockwise by 90 degrees.
rotateCW :: Floating a => Vec2 a -> Vec2 a
rotateCW (Vec2 x y) = Vec2 y (-x)


