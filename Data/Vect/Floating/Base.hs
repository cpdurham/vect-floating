{-# LANGUAGE DeriveDataTypeable, MultiParamTypeClasses, FunctionalDependencies, GeneralizedNewtypeDeriving, StandaloneDeriving, GADTs, FlexibleInstances, FlexibleContexts, ScopedTypeVariables #-}

module Data.Vect.Floating.Base
  ( AbelianGroup(..) , vecSum
  , MultSemiGroup(..) , Ring , semigroupProduct
  , LeftModule(..) , RightModule(..)
  , Vector(..) , DotProd(..) , CrossProd(..)
  , normalize , distance , angle , angle'
  , UnitVector(..)
  , Pointwise(..)
  , Extend(..) , HasCoordinates(..) , Dimension(..)
  , Matrix(..) , Tensor(..) , Diagonal (..) , Determinant(..)
  , Orthogonal(..) , Projective(..) , MatrixNorms(..)
  , Vec2(..) , Vec3(..) , Vec4(..)
  , Mat2(..) , Mat3(..) , Mat4(..)
  , Ortho2 , Ortho3 , Ortho4
  , Normal2 , Normal3 , Normal4
  , Proj3 , Proj4
  , mkVec2 , mkVec3 , mkVec4
  , project , project' , projectUnsafe , flipNormal
  , householder, householderOrtho
  )
  where

import Data.Typeable
import Control.Monad
import System.Random  
import Foreign

--------------------------------------------------------------------------------
-- class declarations

class AbelianGroup g where
  (&+) :: g -> g -> g
  (&-) :: g -> g -> g
  neg  :: g -> g
  zero :: g

infixl 6 &+
infixl 6 &- 

vecSum :: AbelianGroup g => [g] -> g
vecSum l = foldl (&+) zero l 

class MultSemiGroup r where
  (.*.) :: r -> r -> r
  one   :: r

class (AbelianGroup r, MultSemiGroup r) => Ring r 

infixl 7 .*. 

-- was: ringProduct :: Ring r => [r] -> r
semigroupProduct :: MultSemiGroup r => [r] -> r 
semigroupProduct l = foldl (.*.) one l

class LeftModule r m where
  lmul :: r -> m -> m
  (*.) :: r -> m -> m
  (*.) = lmul

class RightModule m r | m -> r, r -> m where
  rmul :: m -> r -> m
  (.*) :: m -> r -> m
  (.*) = rmul

-- I'm not really sure about this.. may actually degrade the performance in some cases?  
{- RULES
"matrix multiplication left"   forall m n x.  (n .*. m) *. x = n *. (m *. x)  
"matrix multiplication right"  forall m n x.  x .* (m .*. n) = (x .* m) .* n
  -}

infixr 7 *.
infixl 7 .*

class AbelianGroup (v a) => Vector a v where
  mapVec    :: (a -> a) -> v a -> v a
  scalarMul :: a -> v a -> v a
  (*&) ::      a -> v a -> v a
  (&*) ::      v a -> a -> v a 
  (*&) s v = scalarMul s v
  (&*) v s = scalarMul s v

infixr 7 *&
infixl 7 &*

{-# RULES
"scalar multiplication left"   forall (s :: Num s => s) (t :: Num t => t) x. t *& (s *& x) = (t*s) *& x 
"scalar multiplication right"  forall (s :: Num s => s) (t :: Num t => t) x.  (x &* s) &* t = x &* (s*t)  
  #-}

class Floating a => DotProd a v where
  (&.) :: v a -> v a -> a
  norm    :: v a -> a
  normsqr :: v a -> a
  len     :: v a -> a
  lensqr  :: v a -> a
  len = norm
  lensqr = normsqr
  dotprod :: v a -> v a -> a
  normsqr v = (v &. v)
  norm = sqrt . lensqr
  dotprod = (&.)

infix 7 &.

{-# RULES
"len/square 1"   forall x.  (len x)*(len x) = lensqr x
"len/square 2"   forall x.  (len x)^2 = lensqr x
"norm/square 1"  forall x.  (norm x)*(norm x) = normsqr x
"norm/square 2"  forall x.  (norm x)^2 = normsqr x
  #-}


normalize :: (Vector a v, DotProd a v) => v a -> v a
normalize v = scalarMul (recip (len v)) v

distance :: (Vector a v, DotProd a v) => v a -> v a -> a
distance x y = norm (x &- y)

-- | the angle between two vectors
angle :: (Vector a v, DotProd a v) => v a -> v a -> a 
angle x y = acos $ (x &. y) / (norm x * norm y)

-- | the angle between two unit vectors
angle' {- ' CPP is sensitive to primes -} :: (Vector a v, UnitVector a v u, DotProd a v) => u a -> u a -> a
angle' x y = acos (fromNormal x &. fromNormal y)

{-# RULES
"normalize is idempotent"  forall x. normalize (normalize x) = normalize x
  #-}

class (Vector a v, DotProd a v) => UnitVector a v u | u -> v, v -> u where
  mkNormal         :: v a -> u a       -- ^ normalizes the input
  toNormalUnsafe   :: v a -> u a       -- ^ does not normalize the input!
  fromNormal       :: u a -> v a
  fromNormalRadius :: a -> u a -> v a
  fromNormalRadius t n = t *& fromNormal n 

-- | Projects the first vector down to the hyperplane orthogonal to the second (unit) vector
project' :: (Vector a v, UnitVector a v u, DotProd a v) => v a -> u a -> v a
project' what dir = projectUnsafe what (fromNormal dir)

-- | Direction (second argument) is assumed to be a /unit/ vector!
projectUnsafe :: (Vector a v, DotProd a v) => v a -> v a -> v a
projectUnsafe what dir = what &- dir &* (what &. dir)

project :: (Vector a v, DotProd a v) => v a -> v a -> v a
project what dir = what &- dir &* ((what &. dir) / (dir &. dir))

-- | Since unit vectors are not a group, we need a separate function.
flipNormal :: UnitVector a v n => n a -> n a
flipNormal = toNormalUnsafe . neg . fromNormal 

-- | Cross product
class CrossProd v where
  crossprod :: v -> v -> v
  (&^)      :: v -> v -> v
  (&^) = crossprod
 
-- | Pointwise multiplication 
class Pointwise v where
  pointwise :: v -> v -> v
  (&!)      :: v -> v -> v
  (&!) = pointwise 

infix 7 &^
infix 7 &!

class HasCoordinates v x | v->x where
  _1 :: v -> x
  _2 :: v -> x
  _3 :: v -> x
  _4 :: v -> x

-- | conversion between vectors (and matrices) of different dimensions
class Extend a u v where
  extendZero :: u a -> v a          -- ^ example: @extendZero (Vec2 5 6) = Vec4 5 6 0 0@
  extendWith :: a -> u a -> v a   -- ^ example: @extendWith 1 (Vec2 5 6) = Vec4 5 6 1 1@
  trim :: v a -> u a                -- ^ example: @trim (Vec4 5 6 7 8) = Vec2 5 6@

-- | makes a diagonal matrix from a vector
class Diagonal s t | t->s where
  diag :: s -> t

class Matrix m where
  transpose :: m -> m 
  inverse :: m -> m
  idmtx :: m

{-# RULES
"transpose is an involution"  forall m. transpose (transpose m) = m
"inverse is an involution"    forall m. inverse (inverse m) = m
  #-}
  
class Matrix (m a) => Orthogonal a m o | m -> o, o -> m where
  fromOrtho     :: o a -> m a
  toOrthoUnsafe :: m a -> o a
  
class (AbelianGroup m, Matrix m) => MatrixNorms a m where
  frobeniusNorm  :: m -> a       -- ^ the frobenius norm (= euclidean norm in the space of matrices)
  matrixDistance :: m -> m -> a  -- ^ euclidean distance in the space of matrices
  operatorNorm   :: m -> a      -- ^ (euclidean) operator norm (not implemented yet)
  matrixDistance m n = frobeniusNorm (n &- m)
  operatorNorm = error "operatorNorm: not implemented yet"
  
-- | Outer product (could be unified with Diagonal?)
class Tensor t v | t->v where
  outer :: v -> v -> t
    
class Determinant a m where
  det :: m -> a

class Dimension a where
  dim :: a -> Int
     
-- | Householder matrix, see <http://en.wikipedia.org/wiki/Householder_transformation>.  
-- In plain words, it is the reflection to the hyperplane orthogonal to the input vector.
householder :: (Vector a v, UnitVector a v u, Matrix (m a), Vector a m, Tensor (m a) (v a)) => u a -> m a
householder u = idmtx &- (2 *& outer v v) 
  where v = fromNormal u

householderOrtho :: (Vector a v, UnitVector a v u, Matrix (m a), Vector a m, Tensor (m a) (v a), Orthogonal a m o) => u a -> o a
householderOrtho = toOrthoUnsafe . householder

-- | \"Projective\" matrices have the following form: the top left corner
-- is an any matrix, the bottom right corner is 1, and the top-right
-- column is zero. These describe the affine orthogonal transformation of
-- the space one dimension less.
class (Vector a v, Orthogonal a n o, Diagonal (v a) (n a)) => Projective a v n o m p | m -> p, p -> m, p -> o, o -> p, p -> n, n -> p, p -> v, v -> p, n -> o, n -> v, v -> n where
  fromProjective     :: p a -> m a
  toProjectiveUnsafe :: m a -> p a
  orthogonal         :: o a -> p a
  linear             :: n a -> p a
  translation        :: v a -> p a
  scaling            :: v a -> p a

--------------------------------------------------------------------------------
-- Vec / Mat datatypes
  
data Vec2 a = Vec2 !a !a 
  deriving (Read,Show,Typeable)
data Vec3 a = Vec3 !a !a !a
  deriving (Read,Show,Typeable)
data Vec4 a = Vec4 !a !a !a !a
  deriving (Read,Show,Typeable)

-- | The components are /row/ vectors 

data Mat2 a = Mat2 !(Vec2 a) !(Vec2 a) deriving (Read,Show)
data Mat3 a = Mat3 !(Vec3 a) !(Vec3 a) !(Vec3 a) deriving (Read,Show)
data Mat4 a = Mat4 !(Vec4 a) !(Vec4 a) !(Vec4 a) !(Vec4 a) deriving (Read,Show)

-- | The assumption when dealing with these is always that they are of unit length.
-- Also, interpolation works differently.
newtype Normal2 a = Normal2 (Vec2 a) deriving (Read,Show,Storable,Dimension,Typeable) 
newtype Normal3 a = Normal3 (Vec3 a) deriving (Read,Show,Storable,Dimension,Typeable)
newtype Normal4 a = Normal4 (Vec4 a) deriving (Read,Show,Storable,Dimension,Typeable)

deriving instance Floating a => DotProd a Normal2
deriving instance Floating a => DotProd a Normal3
deriving instance Floating a => DotProd a Normal4

mkVec2 :: (a,a) -> Vec2 a
mkVec3 :: (a,a,a) -> Vec3 a
mkVec4 :: (a,a,a,a) -> Vec4 a

mkVec2 (x,y)     = Vec2 x y
mkVec3 (x,y,z)   = Vec3 x y z
mkVec4 (x,y,z,w) = Vec4 x y z w

-- | Orthogonal matrices.
--
-- Note: the "Random" instances generates orthogonal matrices with determinant 1
-- (that is, orientation-preserving orthogonal transformations)!
newtype Ortho2 a = Ortho2 (Mat2 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)
newtype Ortho3 a = Ortho3 (Mat3 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)
newtype Ortho4 a = Ortho4 (Mat4 a) deriving (Read,Show,Storable,MultSemiGroup,Determinant a,Dimension)

-- | Projective matrices, encoding affine transformations in dimension one less.
newtype Proj3 a = Proj3 (Mat3 a) deriving (Read,Show,Storable,MultSemiGroup)
newtype Proj4 a = Proj4 (Mat4 a) deriving (Read,Show,Storable,MultSemiGroup)

--------------------------------------------------------------------------------
-- Unit vectors
  
instance Floating a => UnitVector a Vec2 Normal2 where
  mkNormal v = Normal2 (normalize v)
  fromNormal (Normal2 v) = v 
  toNormalUnsafe = Normal2

instance Floating a => UnitVector a Vec3 Normal3 where
  mkNormal v = Normal3 (normalize v)
  fromNormal (Normal3 v) = v 
  toNormalUnsafe = Normal3

instance Floating a => UnitVector a Vec4 Normal4 where
  mkNormal v = Normal4 (normalize v)
  fromNormal (Normal4 v) = v 
  toNormalUnsafe = Normal4

_rndUnit :: (Ord a, RandomGen g, Random (v a), Vector a v, DotProd a v) => g -> (v a,g)
_rndUnit g = 
  if d > 0.01
    then ( v &* (1.0/d) , h )
    else _rndUnit h
  where
    (v,h) = random g
    d = norm v
    
instance (Floating a, Random a, Ord a) => Random (Normal2 a) where
  random g = let (v,h) = _rndUnit g in (Normal2 v, h)  
  randomR _ = random

instance (Floating a, Random a, Ord a) => Random (Normal3 a) where
  random g = let (v,h) = _rndUnit g in (Normal3 v, h)  
  randomR _ = random

instance (Floating a, Random a, Ord a) => Random (Normal4 a) where
  random g = let (v,h) = _rndUnit g in (Normal4 v, h)  
  randomR _ = random

instance Floating a => CrossProd (Normal3 a) where
  crossprod (Normal3 v) (Normal3 w) = mkNormal (crossprod v w)

--------------------------------------------------------------------------------
-- Orthogonal matrices

instance Floating a =>  Orthogonal a Mat2 Ortho2 where
  fromOrtho (Ortho2 o) = o
  toOrthoUnsafe = Ortho2

instance Floating a => Orthogonal a Mat3 Ortho3 where
  fromOrtho (Ortho3 o) = o
  toOrthoUnsafe = Ortho3 

instance Floating a => Orthogonal a Mat4 Ortho4 where
  fromOrtho (Ortho4 o) = o
  toOrthoUnsafe = Ortho4

------

instance Floating a => Matrix (Ortho2 a) where
  transpose (Ortho2 o) = Ortho2 (transpose o)
  idmtx = Ortho2 idmtx
  inverse = transpose

instance Floating a => Matrix (Ortho3 a) where
  transpose (Ortho3 o) = Ortho3 (transpose o)
  idmtx = Ortho3 idmtx
  inverse = transpose

instance Floating a => Matrix (Ortho4 a) where
  transpose (Ortho4 o) = Ortho4 (transpose o)
  idmtx = Ortho4 idmtx
  inverse = transpose

------

instance (Floating a, Ord a, Random a) => Random (Ortho2 a) where
  random g = let (o,h) = _rndOrtho2 g in (toOrthoUnsafe (_flip1stRow2 o), h)
  randomR _ = random

instance (Floating a, Ord a, Random a) => Random (Ortho3 a) where
  random g = let (o,h) = _rndOrtho3 g in (toOrthoUnsafe (             o), h)
  randomR _ = random

instance (Floating a, Ord a, Random a) => Random (Ortho4 a) where
  random g = let (o,h) = _rndOrtho4 g in (toOrthoUnsafe (_flip1stRow4 o), h)
  randomR _ = random

------

-- determinant will be -1
_rndOrtho2 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (Mat2 a, g)
_rndOrtho2 g = (h2, g1) where
  h2 = householder u2
  (u2,g1) = random g   

-- generates a uniformly random orthogonal 3x3 matrix 
-- /with determinant +1/, with respect to the Haar measure of SO3.
--
-- see Theorem 4 in:
-- Francesco Mezzadri: How to Generate Random Matrices from the Classical Compact Groups 
-- Notices of the AMS, May 2007 issue
-- <http://www.ams.org/notices/200705/fea-mezzadri-web.ps>
_rndOrtho3 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (Mat3 a, g) 
_rndOrtho3 g = ( (h3 .*. m3), g2) where
  m3 = (extendWith :: Floating a => a -> Mat2 a -> Mat3 a) 1 o2 
  h3 = householder u3
  (u3,g1) = random g
  (o2,g2) = _rndOrtho2 g1

-- determinant will be -1
_rndOrtho4 :: (Floating a, Random a, Ord a, RandomGen g) => g -> (Mat4 a, g) 
_rndOrtho4 g = ( (h4 .*. m4), g2) where
  m4 = (extendWith :: Floating a => a -> Mat3 a -> Mat4 a) 1 o3 
  h4 = householder u4
  (u4,g1) = random g
  (o3,g2) = _rndOrtho3 g1

------

_flip1stRow2 :: Floating a => Mat2 a -> Mat2 a
_flip1stRow2 (Mat2 a b) = Mat2 (neg a) b

_flip1stRow3 :: Floating a => Mat3 a -> Mat3 a
_flip1stRow3 (Mat3 a b c) = Mat3 (neg a) b c

_flip1stRow4 :: Floating a => Mat4 a -> Mat4 a
_flip1stRow4 (Mat4 a b c d) = Mat4 (neg a) b c d

--------------------------------------------------------------------------------
-- projective matrices
  
instance Floating a => Projective a Vec2 Mat2 Ortho2 Mat3 Proj3 where
  fromProjective (Proj3 m) = m
  toProjectiveUnsafe = Proj3
  orthogonal = Proj3 . extendWith (1 :: a) . fromOrtho
  linear     = Proj3 . extendWith (1 :: a)
  translation v = Proj3 $ Mat3 (Vec3 1 0 0) (Vec3 0 1 0) (extendWith (1 :: a) v)
  scaling     v = Proj3 $ diag (extendWith (1 :: a) v)
  
instance Floating a => Projective a Vec3 Mat3 Ortho3 Mat4 Proj4 where
  fromProjective (Proj4 m) = m
  toProjectiveUnsafe = Proj4
  orthogonal = Proj4 . extendWith (1 :: a) . fromOrtho 
  linear     = Proj4 . extendWith (1 :: a)
  translation v = Proj4 $ Mat4 (Vec4 1 0 0 0) (Vec4 0 1 0 0) (Vec4 0 0 1 0) (extendWith (1 :: a) v)
  scaling     v = Proj4 $ diag (extendWith (1 :: a) v)

instance Floating a => Matrix (Proj3 a) where
  idmtx = Proj3 idmtx
  transpose (Proj3 m) = Proj3 (transpose m)
  inverse = _invertProj3

instance Floating a => Matrix (Proj4 a) where
  idmtx = Proj4 idmtx
  transpose (Proj4 m) = Proj4 (transpose m)
  inverse = _invertProj4

_invertProj3 :: (Floating a, Extend a Vec2 Vec3) => Proj3 a -> Proj3 a
_invertProj3 (Proj3 mat@(Mat3 _ _ t)) = 
  Proj3 $ Mat3 (extendZero a) (extendZero b) (extendWith 1 t')
  where
    t' = neg $ (trim :: Extend a Vec2 Vec3 => Vec3 a -> Vec2 a) t .* invm2
    invm2@(Mat2 a b) = inverse $ trim mat

-- Inverts a projective 4x4 matrix. But you can simply use "inverse" instead.
-- We assume that the bottom-right corner is 1.
_invertProj4 :: Floating a => Proj4 a -> Proj4 a
_invertProj4 (Proj4 mat@(Mat4 _ _ _ t)) = 
  Proj4 $ Mat4 (extendZero a) (extendZero b) (extendZero c) (extendWith 1 t') 
  where
    t' = neg $ (trim :: Extend a Vec3 Vec4 => Vec4 a -> Vec3 a) t .* invm3 
    invm3@(Mat3 a b c) = inverse $ trim mat

--------------------------------------------------------------------------------
-- Vec2 instances

instance Floating a => HasCoordinates (Vec2 a) a where
  _1 (Vec2 x _) = x
  _2 (Vec2 _ y) = y
  _3 _ = error "has only 2 coordinates"
  _4 _ = error "has only 2 coordinates"

instance Floating a => AbelianGroup (Vec2 a) where
  (&+) (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1+x2) (y1+y2) 
  (&-) (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1-x2) (y1-y2)
  neg  (Vec2 x y)                = Vec2 (-x) (-y)
  zero = Vec2 0 0
  
instance Floating a => Vector a Vec2 where
  scalarMul s (Vec2 x y) = Vec2 (s*x) (s*y)
  mapVec    f (Vec2 x y) = Vec2 (f x) (f y)
  
instance Floating a => DotProd a Vec2 where
  (&.) (Vec2 x1 y1) (Vec2 x2 y2) = x1*x2 + y1*y2

instance Floating a => Pointwise (Vec2 a) where
  pointwise (Vec2 x1 y1) (Vec2 x2 y2) = Vec2 (x1*x2) (y1*y2)

instance Floating a => Determinant a (Vec2 a, Vec2 a) where
  det (Vec2 x1 y1 , Vec2 x2 y2) = x1*y2 - x2*y1

{-     
instance Show Vec2 where
  show (Vec2 x y) = "( " ++ show x ++ " , " ++ show y ++ " )"
-}

instance (Floating a, Random a) => Random (Vec2 a) where
  random = randomR (Vec2 (-1) (-1),Vec2 1 1)
  randomR (Vec2 a b, Vec2 c d) gen = 
    let (x,gen1) = randomR (a,c) gen
        (y,gen2) = randomR (b,d) gen1
    in (Vec2 x y, gen2)
     
instance (Floating a, Storable a) => Storable (Vec2 a) where
  sizeOf    _ = 2 * sizeOf (undefined :: a)
  alignment _ = sizeOf (undefined :: a)
  
  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p 
    y <- peekByteOff p k
    return (Vec2 x y)
    
  poke q (Vec2 x y) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p   x
    pokeByteOff p k y

instance Floating a => Dimension (Vec2 a) where dim _ = 2

--------------------------------------------------------------------------------                    
-- Mat2 instances

instance Floating a => HasCoordinates (Mat2 a) (Vec2 a) where
  _1 (Mat2 x _) = x
  _2 (Mat2 _ y) = y
  _3 _ = error "has only 2 coordinates"
  _4 _ = error "has only 2 coordinates"

instance Floating a => Matrix (Mat2 a) where
  transpose (Mat2 row1 row2) = 
    Mat2 (Vec2 (_1 row1) (_1 row2)) 
         (Vec2 (_2 row1) (_2 row2)) 
  idmtx = Mat2 (Vec2 1 0) (Vec2 0 1)
  inverse (Mat2 (Vec2 a b) (Vec2 c d)) = 
    Mat2 (Vec2 (d*r) (-b*r)) (Vec2 (-c*r) (a*r)) 
    where r = 1.0 / (a*d - b*c)

instance Floating a => AbelianGroup (Mat2 a) where
  (&+) (Mat2 r1 r2) (Mat2 s1 s2) = Mat2 (r1 &+ s1) (r2 &+ s2)
  (&-) (Mat2 r1 r2) (Mat2 s1 s2) = Mat2 (r1 &- s1) (r2 &- s2)
  neg  (Mat2 r1 r2)              = Mat2 (neg r1) (neg r2)  
  zero = Mat2 zero zero  
  
instance Floating a => Vector a Mat2 where
  scalarMul s (Mat2 r1 r2) = Mat2 (g r1) (g r2) where g = scalarMul s
  mapVec    f (Mat2 r1 r2) = Mat2 (g r1) (g r2) where g = mapVec f

instance Floating a => MultSemiGroup (Mat2 a) where
  (.*.) (Mat2 r1 r2) n = 
    let (Mat2 c1 c2) = transpose n
    in Mat2 (Vec2 (r1 &. c1) (r1 &. c2))
            (Vec2 (r2 &. c1) (r2 &. c2))
  one = idmtx 

instance Floating a => Ring (Mat2 a)

instance Floating a => LeftModule (Mat2 a) (Vec2 a) where
  lmul (Mat2 row1 row2) v = Vec2 (row1 &. v) (row2 &. v) 
  
instance Floating a => RightModule (Vec2 a) (Mat2 a) where
  rmul v mt = lmul (transpose mt) v

instance Floating a => Diagonal (Vec2 a) (Mat2 a) where
  diag (Vec2 x y) = Mat2 (Vec2 x 0) (Vec2 0 y)

instance Floating a => Tensor (Mat2 a) (Vec2 a) where
  outer (Vec2 a b) (Vec2 x y) = Mat2
    (Vec2 (a*x) (a*y))
    (Vec2 (b*x) (b*y))

instance Floating a => Determinant a (Mat2 a) where
  det (Mat2 (Vec2 a b) (Vec2 c d)) = a*d - b*c 

{-
instance Show Mat2 where
  show (Mat2 r1 r2) = show r1 ++ "\n" ++ show r2
-}

instance (Floating a, Storable a) => Storable (Mat2 a) where
  sizeOf    _ = 2 * sizeOf (undefined :: Vec2 a)
  alignment _ = alignment  (undefined :: Vec2 a)
  
  peek q = do
    let p = castPtr q :: Ptr (Vec2 a)
        k = sizeOf (undefined :: Vec2 a)
    r1 <- peek        p 
    r2 <- peekByteOff p k
    return (Mat2 r1 r2)
    
  poke q (Mat2 r1 r2) = do
    let p = castPtr q :: Ptr (Vec2 a)
        k = sizeOf (undefined :: Vec2 a)
    poke        p   r1
    pokeByteOff p k r2

instance (Floating a, Random a) => Random (Mat2 a) where
  random = randomR (Mat2 v1 v1 , Mat2 v2 v2) where 
    v1 = Vec2 (-1) (-1) 
    v2 = Vec2   1    1
  randomR (Mat2 a b, Mat2 c d) gen = 
    let (x,gen1) = randomR (a,c) gen
        (y,gen2) = randomR (b,d) gen1
    in (Mat2 x y, gen2)
          
instance Floating a => Dimension (Mat2 a) where dim _ = 2
     
instance (Floating a) => MatrixNorms a (Mat2 a) where 
  frobeniusNorm (Mat2 r1 r2) =  
    sqrt $
      normsqr r1 + 
      normsqr r2
     
instance Floating a => Pointwise (Mat2 a) where
  pointwise (Mat2 x1 y1) (Mat2 x2 y2) = Mat2 (x1 &! x2) (y1 &! y2)
       
--------------------------------------------------------------------------------     
-- Vec3 instances

instance Floating a => HasCoordinates (Vec3 a) a where
  _1 (Vec3 x _ _) = x
  _2 (Vec3 _ y _) = y
  _3 (Vec3 _ _ z) = z
  _4 _ = error "has only 3 coordinates"

instance Floating a => AbelianGroup (Vec3 a) where
  (&+) (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1+x2) (y1+y2) (z1+z2) 
  (&-) (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1-x2) (y1-y2) (z1-z2) 
  neg  (Vec3 x y z)                    = Vec3 (-x) (-y) (-z)
  zero = Vec3 0 0 0
  
instance Floating a => Vector a Vec3 where
  scalarMul s (Vec3 x y z) = Vec3 (s*x) (s*y) (s*z)
  mapVec    f (Vec3 x y z) = Vec3 (f x) (f y) (f z)

instance Floating a => DotProd a Vec3 where
  (&.) (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = x1*x2 + y1*y2 + z1*z2

instance Num a => Pointwise (Vec3 a) where
  pointwise (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (x1*x2) (y1*y2) (z1*z2)

{-
instance Show Vec3 where
  show (Vec3 x y z) = "( " ++ show x ++ " , " ++ show y ++ " , " ++ show z ++ " )"
-}

instance (Floating a, Random a) => Random (Vec3 a) where
  random = randomR (Vec3 (-1) (-1) (-1),Vec3 1 1 1)
  randomR (Vec3 a b c, Vec3 d e f) gen = 
    let (x,gen1) = randomR (a,d) gen
        (y,gen2) = randomR (b,e) gen1
        (z,gen3) = randomR (c,f) gen2  
    in (Vec3 x y z, gen3)
      
instance Floating a => CrossProd (Vec3 a) where
  crossprod (Vec3 x1 y1 z1) (Vec3 x2 y2 z2) = Vec3 (y1*z2-y2*z1) (z1*x2-z2*x1) (x1*y2-x2*y1) 

instance Floating a => Determinant a (Vec3 a, Vec3 a, Vec3 a) where
  det (u,v,w) = u &. (v &^ w)  
 
instance (Floating a, Storable a) => Storable (Vec3 a) where
  sizeOf    _ = 3 * sizeOf (undefined :: a)
  alignment _ = sizeOf (undefined :: a)
  
  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p 
    y <- peekByteOff p (k  )
    z <- peekByteOff p (k+k)
    return (Vec3 x y z)
    
  poke q (Vec3 x y z) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       x
    pokeByteOff p (k  ) y
    pokeByteOff p (k+k) z

instance Floating a => Dimension (Vec3 a) where dim _ = 3

--------------------------------------------------------------------------------   
-- Mat3 instances

instance Floating a => HasCoordinates (Mat3 a) (Vec3 a) where
  _1 (Mat3 x _ _) = x
  _2 (Mat3 _ y _) = y
  _3 (Mat3 _ _ z) = z
  _4 _ = error "has only 3 coordinates"  

instance Floating a => Matrix (Mat3 a) where

  transpose (Mat3 row1 row2 row3) = 
    Mat3 (Vec3 (_1 row1) (_1 row2) (_1 row3)) 
         (Vec3 (_2 row1) (_2 row2) (_2 row3)) 
         (Vec3 (_3 row1) (_3 row2) (_3 row3)) 
         
  idmtx = Mat3 (Vec3 1 0 0) (Vec3 0 1 0) (Vec3 0 0 1)
  
  inverse (Mat3 (Vec3 a b c) (Vec3 e f g) (Vec3 i j k)) = 
    Mat3 (Vec3 (d11*r) (d21*r) (d31*r))  
         (Vec3 (d12*r) (d22*r) (d32*r))  
         (Vec3 (d13*r) (d23*r) (d33*r))  
    where
      r = 1.0 / ( a*d11 + b*d12 + c*d13 )

      d11 = f*k - g*j
      d12 = g*i - e*k
      d13 = e*j - f*i

      d31 = b*g - c*f
      d32 = c*e - a*g
      d33 = a*f - b*e

      d21 = c*j - b*k 
      d22 = a*k - c*i 
      d23 = b*i - a*j 

instance Floating a => AbelianGroup (Mat3 a) where
  (&+) (Mat3 r1 r2 r3) (Mat3 s1 s2 s3) = Mat3 (r1 &+ s1) (r2 &+ s2) (r3 &+ s3)
  (&-) (Mat3 r1 r2 r3) (Mat3 s1 s2 s3) = Mat3 (r1 &- s1) (r2 &- s2) (r3 &- s3)
  neg  (Mat3 r1 r2 r3)                 = Mat3 (neg r1) (neg r2) (neg r3) 
  zero = Mat3 zero zero zero 

instance Floating a => Vector a Mat3 where
  scalarMul s (Mat3 r1 r2 r3) = Mat3 (g r1) (g r2) (g r3) where g = scalarMul s
  mapVec    f (Mat3 r1 r2 r3) = Mat3 (g r1) (g r2) (g r3) where g = mapVec f

instance Floating a => MultSemiGroup (Mat3 a) where
  (.*.) (Mat3 r1 r2 r3) n = 
    let (Mat3 c1 c2 c3) = transpose n
    in Mat3 (Vec3 (r1 &. c1) (r1 &. c2) (r1 &. c3))
            (Vec3 (r2 &. c1) (r2 &. c2) (r2 &. c3))
            (Vec3 (r3 &. c1) (r3 &. c2) (r3 &. c3))
  one = idmtx 

instance Floating a => Ring (Mat3 a)

instance Floating a => LeftModule (Mat3 a) (Vec3 a) where
  lmul (Mat3 row1 row2 row3) v = Vec3 (row1 &. v) (row2 &. v) (row3 &. v)
  
instance Floating a => RightModule (Vec3 a) (Mat3 a) where
  rmul v mt = lmul (transpose mt) v

instance Floating a => Diagonal (Vec3 a) (Mat3 a) where
  diag (Vec3 x y z) = Mat3 (Vec3 x 0 0) (Vec3 0 y 0) (Vec3 0 0 z)

instance Floating a => Tensor (Mat3 a) (Vec3 a) where
  outer (Vec3 a b c) (Vec3 x y z) = Mat3
    (Vec3 (a*x) (a*y) (a*z))
    (Vec3 (b*x) (b*y) (b*z))
    (Vec3 (c*x) (c*y) (c*z))

instance Floating a => Determinant a (Mat3 a) where
  det (Mat3 r1 r2 r3) = det (r1,r2,r3)

{-
instance Show Mat3 where
  show (Mat3 r1 r2 r3) = show r1 ++ "\n" ++ show r2 ++ "\n" ++ show r3
-}

instance (Floating a, Storable a) => Storable (Mat3 a) where
  sizeOf    _ = 3 * sizeOf (undefined::Vec3 a)
  alignment _ = alignment  (undefined::Vec3 a)
  
  peek q = do
    let p = castPtr q :: Ptr (Vec3 a)
        k = sizeOf (undefined::Vec3 a)
    r1 <- peek        p 
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    return (Mat3 r1 r2 r3)
    
  poke q (Mat3 r1 r2 r3) = do
    let p = castPtr q :: Ptr (Vec3 a)
        k = sizeOf (undefined::Vec3 a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3

instance (Floating a, Random a) => Random (Mat3 a) where
  random = randomR (Mat3 v1 v1 v1 , Mat3 v2 v2 v2) where
    v1 = Vec3 (-1) (-1) (-1)
    v2 = Vec3   1    1    1
  randomR (Mat3 a b c, Mat3 d e f) gen = 
    let (x,gen1) = randomR (a,d) gen
        (y,gen2) = randomR (b,e) gen1
        (z,gen3) = randomR (c,f) gen2  
    in (Mat3 x y z, gen3)
   
instance Floating a => Dimension (Mat3 a) where dim _ = 3
  
instance Floating a => MatrixNorms a (Mat3 a) where 
  frobeniusNorm (Mat3 r1 r2 r3)  = 
    sqrt $
      normsqr r1 + 
      normsqr r2 + 
      normsqr r3 

instance Floating a => Pointwise (Mat3 a) where
  pointwise (Mat3 x1 y1 z1) (Mat3 x2 y2 z2) = Mat3 (x1 &! x2) (y1 &! y2) (z1 &! z2)
    
--------------------------------------------------------------------------------
-- Vec4 instances

instance Floating a => HasCoordinates (Vec4 a) a where
  _1 (Vec4 x _ _ _) = x
  _2 (Vec4 _ y _ _) = y
  _3 (Vec4 _ _ z _) = z
  _4 (Vec4 _ _ _ w) = w

instance Floating a => AbelianGroup (Vec4 a) where
  (&+) (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = Vec4 (x1+x2) (y1+y2) (z1+z2) (w1+w2)
  (&-) (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = Vec4 (x1-x2) (y1-y2) (z1-z2) (w1-w2)
  neg  (Vec4 x y z w)                        = Vec4 (-x) (-y) (-z) (-w)
  zero = Vec4 0 0 0 0
  
instance Floating a => Vector a Vec4 where
  scalarMul s (Vec4 x y z w) = Vec4 (s*x) (s*y) (s*z) (s*w)
  mapVec    f (Vec4 x y z w) = Vec4 (f x) (f y) (f z) (f w)

instance Floating a => DotProd a Vec4 where
  (&.) (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = x1*x2 + y1*y2 + z1*z2 + w1*w2

instance Num a => Pointwise (Vec4 a) where
  pointwise (Vec4 x1 y1 z1 w1) (Vec4 x2 y2 z2 w2) = Vec4 (x1*x2) (y1*y2) (z1*z2) (w1*w2)

{-
instance Show Vec4 where
  show (Vec4 x y z w) = "( " ++ show x ++ " , " ++ show y ++ " , " ++ show z ++ " , " ++ show w ++ " )"
-}

instance (Floating a, Random a) => Random (Vec4 a) where
  random = randomR (Vec4 (-1) (-1) (-1) (-1),Vec4 1 1 1 1)
  randomR (Vec4 a b c d, Vec4 e f g h) gen = 
    let (x,gen1) = randomR (a,e) gen
        (y,gen2) = randomR (b,f) gen1
        (z,gen3) = randomR (c,g) gen2  
        (w,gen4) = randomR (d,h) gen3  
    in (Vec4 x y z w, gen4)
           
instance (Floating a, Storable a) => Storable (Vec4 a) where
  sizeOf    _ = 4 * sizeOf (undefined :: a)
  alignment _ = sizeOf (undefined :: a)
  
  peek q = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    x <- peek        p 
    y <- peekByteOff p (k  )
    z <- peekByteOff p (k+k)
    w <- peekByteOff p (3*k)
    return (Vec4 x y z w)
    
  poke q (Vec4 x y z w) = do
    let p = castPtr q :: Ptr a
        k = sizeOf (undefined :: a)
    poke        p       x
    pokeByteOff p (k  ) y
    pokeByteOff p (k+k) z
    pokeByteOff p (3*k) w

instance Floating a => Dimension (Vec4 a) where dim _ = 4

--------------------------------------------------------------------------------
-- Mat4 instances

instance Floating a => HasCoordinates (Mat4 a) (Vec4 a) where
  _1 (Mat4 x _ _ _) = x
  _2 (Mat4 _ y _ _) = y
  _3 (Mat4 _ _ z _) = z
  _4 (Mat4 _ _ _ w) = w

instance Floating a => Matrix (Mat4 a) where
  transpose (Mat4 row1 row2 row3 row4) = 
    Mat4 (Vec4 (_1 row1) (_1 row2) (_1 row3) (_1 row4)) 
         (Vec4 (_2 row1) (_2 row2) (_2 row3) (_2 row4)) 
         (Vec4 (_3 row1) (_3 row2) (_3 row3) (_3 row4)) 
         (Vec4 (_4 row1) (_4 row2) (_4 row3) (_4 row4)) 
  idmtx = Mat4 (Vec4 1 0 0 0) (Vec4 0 1 0 0) (Vec4 0 0 1 0) (Vec4 0 0 0 1)
  inverse = error "inverse/Mat4: not implemented yet"

instance Floating a => AbelianGroup (Mat4 a) where
  (&+) (Mat4 r1 r2 r3 r4) (Mat4 s1 s2 s3 s4) = Mat4 (r1 &+ s1) (r2 &+ s2) (r3 &+ s3) (r4 &+ s4)
  (&-) (Mat4 r1 r2 r3 r4) (Mat4 s1 s2 s3 s4) = Mat4 (r1 &- s1) (r2 &- s2) (r3 &- s3) (r4 &- s4)
  neg  (Mat4 r1 r2 r3 r4)                    = Mat4 (neg r1) (neg r2) (neg r3) (neg r4) 
  zero = Mat4 zero zero zero zero
  
instance Floating a => Vector a Mat4 where
  scalarMul s (Mat4 r1 r2 r3 r4) = Mat4 (g r1) (g r2) (g r3) (g r4) where g = scalarMul s
  mapVec    f (Mat4 r1 r2 r3 r4) = Mat4 (g r1) (g r2) (g r3) (g r4) where g = mapVec f

instance Floating a => MultSemiGroup (Mat4 a) where
  (.*.) (Mat4 r1 r2 r3 r4) n = 
    let (Mat4 c1 c2 c3 c4) = transpose n
    in Mat4 (Vec4 (r1 &. c1) (r1 &. c2) (r1 &. c3) (r1 &. c4))
            (Vec4 (r2 &. c1) (r2 &. c2) (r2 &. c3) (r2 &. c4))
            (Vec4 (r3 &. c1) (r3 &. c2) (r3 &. c3) (r3 &. c4))
            (Vec4 (r4 &. c1) (r4 &. c2) (r4 &. c3) (r4 &. c4))
  one = idmtx 

instance Floating a => Ring (Mat4 a)

instance Floating a => LeftModule (Mat4 a) (Vec4 a) where
  lmul (Mat4 row1 row2 row3 row4) v = Vec4 (row1 &. v) (row2 &. v) (row3 &. v) (row4 &. v)
  
instance Floating a => RightModule (Vec4 a) (Mat4 a) where
  rmul v mt = lmul (transpose mt) v

instance Floating a => Diagonal (Vec4 a) (Mat4 a) where
  diag (Vec4 x y z w) = Mat4 (Vec4 x 0 0 0) (Vec4 0 y 0 0) (Vec4 0 0 z 0) (Vec4 0 0 0 w)

instance Floating a => Tensor (Mat4 a) (Vec4 a) where
  outer (Vec4 a b c d) (Vec4 x y z w) = Mat4
    (Vec4 (a*x) (a*y) (a*z) (a*w))
    (Vec4 (b*x) (b*y) (b*z) (b*w))
    (Vec4 (c*x) (c*y) (c*z) (c*w))
    (Vec4 (d*x) (d*y) (d*z) (d*w))

instance Floating a => Determinant a (Mat4 a) where
  det = error "det/Mat4: not implemented yet" 
  -- det (Mat4 r1 r2 r3 r4) = 

{-
instance Show Mat4 where
  show (Mat4 r1 r2 r3 r4) = show r1 ++ "\n" ++ show r2 ++ "\n" ++ show r3 ++ "\n" ++ show r4
-}

instance (Floating a, Storable a) => Storable (Mat4 a) where
  sizeOf    _ = 4 * sizeOf (undefined::Vec4 a)
  alignment _ = alignment  (undefined::Vec4 a)
  
  peek q = do
    let p = castPtr q :: Ptr (Vec4 a)
        k = sizeOf (undefined :: Vec4 a)
    r1 <- peek        p 
    r2 <- peekByteOff p (k  )
    r3 <- peekByteOff p (k+k)
    r4 <- peekByteOff p (3*k)
    return (Mat4 r1 r2 r3 r4)
    
  poke q (Mat4 r1 r2 r3 r4) = do
    let p = castPtr q :: Ptr (Vec4 a)
        k = sizeOf (undefined :: Vec4 a)
    poke        p       r1
    pokeByteOff p (k  ) r2
    pokeByteOff p (k+k) r3
    pokeByteOff p (3*k) r4

instance (Floating a, Random a) => Random (Mat4 a) where
  random = randomR (Mat4 v1 v1 v1 v1, Mat4 v2 v2 v2 v2) where
    v1 = Vec4 (-1) (-1) (-1) (-1)
    v2 = Vec4   1    1    1    1
  randomR (Mat4 a b c d, Mat4 e f g h) gen = 
    let (x,gen1) = randomR (a,e) gen
        (y,gen2) = randomR (b,f) gen1
        (z,gen3) = randomR (c,g) gen2  
        (w,gen4) = randomR (d,h) gen3  
    in (Mat4 x y z w, gen4)
    
instance Floating a => Dimension (Mat4 a) where dim _ = 4
   
instance Floating a => MatrixNorms a (Mat4 a) where 
  frobeniusNorm (Mat4 r1 r2 r3 r4) = 
    sqrt $
      normsqr r1 + 
      normsqr r2 + 
      normsqr r3 + 
      normsqr r4  
    
instance Floating a => Pointwise (Mat4 a) where
  pointwise (Mat4 x1 y1 z1 w1) (Mat4 x2 y2 z2 w2) = Mat4 (x1 &! x2) (y1 &! y2) (z1 &! z2) (w1 &! w2)
    
--------------------------------------------------------------------------------
-- Extend instances

instance Floating a => Extend a Vec2 Vec3 where
  extendZero   (Vec2 x y) = Vec3 x y 0
  extendWith t (Vec2 x y) = Vec3 x y t
  trim (Vec3 x y _)       = Vec2 x y

instance Floating a => Extend a Vec2 Vec4 where
  extendZero   (Vec2 x y) = Vec4 x y 0 0
  extendWith t (Vec2 x y) = Vec4 x y t t
  trim (Vec4 x y _ _)     = Vec2 x y 

instance Floating a => Extend a Vec3 Vec4 where
  extendZero   (Vec3 x y z) = Vec4 x y z 0
  extendWith t (Vec3 x y z) = Vec4 x y z t
  trim (Vec4 x y z _)       = Vec3 x y z

instance Floating a => Extend a Mat2 Mat3 where
  extendZero   (Mat2 p q) = Mat3 (extendZero p) (extendZero q) zero
  extendWith w (Mat2 p q) = Mat3 (extendZero p) (extendZero q) (Vec3 0 0 w)
  trim (Mat3 p q _) = Mat2 (trim p) (trim q)

instance Floating a => Extend a Mat2 Mat4 where
  extendZero   (Mat2 p q) = Mat4 (extendZero p) (extendZero q) zero zero
  extendWith w (Mat2 p q) = Mat4 (extendZero p) (extendZero q) (Vec4 0 0 w 0) (Vec4 0 0 0 w)
  trim (Mat4 p q _ _) = Mat2 (trim p) (trim q)

instance Floating a => Extend a Mat3 Mat4 where
  extendZero   (Mat3 p q r) = Mat4 (extendZero p) (extendZero q) (extendZero r) zero
  extendWith w (Mat3 p q r) = Mat4 (extendZero p) (extendZero q) (extendZero r) (Vec4 0 0 0 w)
  trim (Mat4 p q r _) = Mat3 (trim p) (trim q) (trim r)
  
--------------------------------------------------------------------------------

  
