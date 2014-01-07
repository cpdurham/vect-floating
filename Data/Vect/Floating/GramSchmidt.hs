{-# LANGUAGE FlexibleInstances #-}

-- | Gram-Schmidt orthogonalization.
-- This module is not re-exported by "Data.Vect".

module Data.Vect.Floating.GramSchmidt 
  ( GramSchmidt(..)
  )
  where

import Data.Vect.Floating.Base

--------------------------------------------------------------------------------

liftPair :: (a -> b) -> (a,a) -> (b,b)
liftPair f (x,y) = (f x, f y)

liftTriple :: (a -> b) -> (a,a,a) -> (b,b,b)
liftTriple f (x,y,z) = (f x, f y, f z)

liftQuadruple :: (a -> b) -> (a,a,a,a) -> (b,b,b,b)
liftQuadruple f (x,y,z,w) = (f x, f y, f z, f w)

--------------------------------------------------------------------------------
    
-- | produces orthogonal\/orthonormal vectors from a set of vectors    
class GramSchmidt a where
  gramSchmidt          :: a -> a   -- ^ does not normalize the vectors!
  gramSchmidtNormalize :: a -> a   -- ^ normalizes the vectors.

{-# RULES
"gramSchmidt is idempotent"  forall a. gramSchmidt (gramSchmidt a) = gramSchmidt a 
"gramSchmidtNormalize is idempotent"  forall a. gramSchmidtNormalize (gramSchmidtNormalize a) = gramSchmidtNormalize a 
  #-}

--------------------------------------------------------------------------------

instance Floating a => GramSchmidt (Vec2 a,Vec2 a) where
  gramSchmidt = gramSchmidtPair
  gramSchmidtNormalize = gramSchmidtNormalizePair
  
instance Floating a => GramSchmidt (Vec3 a,Vec3 a) where
  gramSchmidt = gramSchmidtPair
  gramSchmidtNormalize = gramSchmidtNormalizePair
  
instance Floating a => GramSchmidt (Vec4 a,Vec4 a) where
  gramSchmidt = gramSchmidtPair
  gramSchmidtNormalize = gramSchmidtNormalizePair

----------

instance Floating a => GramSchmidt (Normal2 a,Normal2 a) where
  gramSchmidt          = error "use 'gramSchmidtNormalize' for Normal2!"
  gramSchmidtNormalize = liftPair toNormalUnsafe . gramSchmidtNormalizePair . liftPair fromNormal

instance Floating a => GramSchmidt (Normal3 a,Normal3 a) where
  gramSchmidt          = error "use 'gramSchmidtNormalize' for Normal3!"
  gramSchmidtNormalize = liftPair toNormalUnsafe . gramSchmidtNormalizePair . liftPair fromNormal

instance Floating a => GramSchmidt (Normal4 a,Normal4 a) where
  gramSchmidt          = error "use 'gramSchmidtNormalize' for Normal4!"
  gramSchmidtNormalize = liftPair toNormalUnsafe . gramSchmidtNormalizePair . liftPair fromNormal

----------
  
gramSchmidtPair :: (Vector a v, DotProd a v) => (v a,v a) -> (v a,v a)
gramSchmidtPair (u,v) = (u',v') where 
  u' = u
  v' = project v u'     
  
gramSchmidtNormalizePair :: (Vector a v, DotProd a v) => (v a,v a) -> (v a,v a)
gramSchmidtNormalizePair (u,v) = (u',v') where
  u' = normalize u 
  v' = normalize $ projectUnsafe v u'     

----------

instance Floating a => GramSchmidt (Vec3 a,Vec3 a,Vec3 a) where
  gramSchmidt = gramSchmidtTriple
  gramSchmidtNormalize = gramSchmidtNormalizeTriple
     
instance Floating a => GramSchmidt (Vec4 a,Vec4 a,Vec4 a) where
  gramSchmidt = gramSchmidtTriple
  gramSchmidtNormalize = gramSchmidtNormalizeTriple

instance Floating a => GramSchmidt (Normal3 a,Normal3 a,Normal3 a) where
  gramSchmidt          = error "use 'gramSchmidtNormalize' for Normal3!"
  gramSchmidtNormalize = liftTriple toNormalUnsafe . gramSchmidtNormalizeTriple . liftTriple fromNormal

instance Floating a => GramSchmidt (Normal4 a,Normal4 a,Normal4 a) where
  gramSchmidt          = error "use 'gramSchmidtNormalize' for Normal4!"
  gramSchmidtNormalize = liftTriple toNormalUnsafe . gramSchmidtNormalizeTriple . liftTriple fromNormal

----------

gramSchmidtTriple :: (Vector a v, DotProd a v) => (v a,v a,v a) -> (v a,v a,v a)
gramSchmidtTriple (u,v,w) = (u',v',w') where 
  u' = u
  v' = project v u'     
  w' = project (project w u') v' 
  
gramSchmidtNormalizeTriple :: (Vector a v, DotProd a v) => (v a,v a,v a) -> (v a,v a,v a)
gramSchmidtNormalizeTriple (u,v,w) = (u',v',w') where
  u' = normalize $ u 
  v' = normalize $ projectUnsafe v u'     
  w' = normalize $ projectUnsafe (projectUnsafe w u') v'     

----------

instance Floating a => GramSchmidt (Vec4 a,Vec4 a,Vec4 a,Vec4 a) where
  gramSchmidt          = gramSchmidtQuadruple
  gramSchmidtNormalize = gramSchmidtNormalizeQuadruple 

instance Floating a => GramSchmidt (Normal4 a,Normal4 a,Normal4 a,Normal4 a) where
  gramSchmidt          = error "use 'gramSchmidtNormalize' for Normal4!"
  gramSchmidtNormalize = liftQuadruple toNormalUnsafe . gramSchmidtNormalizeQuadruple . liftQuadruple fromNormal

----------
  
gramSchmidtQuadruple :: (Vector a v, DotProd a v) => (v a,v a,v a,v a) -> (v a,v a,v a,v a)
gramSchmidtQuadruple (u,v,w,z) = (u',v',w',z') where 
  u' = u
  v' = project v u'     
  w' = project (project w u') v' 
  z' = project (project (project z u') v') w'

gramSchmidtNormalizeQuadruple :: (Vector a v, DotProd a v) => (v a,v a,v a,v a) -> (v a,v a,v a,v a)
gramSchmidtNormalizeQuadruple (u,v,w,z) = (u',v',w',z') where
  u' = normalize $ u
  v' = normalize $ projectUnsafe v u'     
  w' = normalize $ projectUnsafe (projectUnsafe w u') v' 
  z' = normalize $ projectUnsafe (projectUnsafe (projectUnsafe z u') v') w'
  
----------
  
