// Learn more about F# at http://fsharp.net
// See the 'F# Tutorial' project for more help.
open System
open System.Numerics
let N= BigInteger 1e9
let rec calcDistance n k a b flag= match flag with 
                                   |true->( (b-1I)*(b-1I)-(a-1I)*(a-1I) ) % N
                                   |false->calcDistance n k (k-a+1I) (k-b+1I) true
let calcCycleDistance n k a b c = (calcDistance n k a c true) + calcDistance n k c a false
let Expectation n k a b c transitionMatrix (distanceVector:BigInteger[,]) = 
    match n with 
    |x when x=1I->(calcDistance n k b a false)+(calcDistance n k a c true)
    |x ->let cycleDistance=((-1I+divTest.fastExpGeneric 2I (n-1I) (fun x y->(x*y) %  N ) ) * (calcCycleDistance n k a b c) )|>divTest.BigIntModK N
         let initDistance= (calcDistance n k b a false)%N
         let fullTransitionMatrix=divTest.fastMatrixExpModK N transitionMatrix (x-1I ) 
         let singletonDistances=( [[1I;0I;0I;0I;0I;0I]]|>array2D|>Matrix.Generic.ofArray2D )*
                                  (fullTransitionMatrix ) *
                                    ( distanceVector|>Matrix.Generic.ofArray2D|>Matrix.Generic.transpose  ) 
         divTest.BigIntModK N (cycleDistance+initDistance+(singletonDistances.[0,0] % N))
                                                                                     
[<EntryPoint>]
let main argv = 
    
    // order (a,b,c) ;(a,c,b);(b,a,c);(b,c,a);(c,b,a);(c,a,b)
    let (transitionMatrix:Matrix<BigInteger>)=Matrix.Generic.ofArray2D (array2D [[0I;1I;1I;0I;0I;0I];[1I;0I;0I;0I;0I;1I];[1I;0I;0I;1I;0I;0I];[0I;0I;1I;0I;1I;0I];[0I;0I;0I;1I;0I;1I];[0I;1I;0I;0I;1I;0I]])
    let calcDistanceVector n k a b c=array2D[[calcDistance n k a c true;calcDistance n k a b true;calcDistance n k b c true ;calcDistance n k b a false; calcDistance n k c a false;calcDistance n k c b false]]
    printfn "%A" (Expectation 3I 1000I 27I 216I 729I transitionMatrix (calcDistanceVector 3I 1000I 27I 216I 729I) )
    let generateA n=divTest.fastExpGeneric 3I n (divTest.BigIntMultModK N )
    let generateB n =divTest.fastExpGeneric 6I n (divTest.BigIntMultModK N )
    let generateC n=divTest.fastExpGeneric 9I n (divTest.BigIntMultModK N )
    let generateK n=divTest.fastExpGeneric 10I n (divTest.BigIntMultModK N )
    let finalAns=List.map (fun n->
                            let k=(generateK n)
                            let a=(generateA n)
                            let b=(generateB n)
                            let c=(generateC n)
                            let currTerm=Expectation n k a b c  transitionMatrix (calcDistanceVector n k a b c) 
                            currTerm) [1I..10000I ]|>List.reduce(fun x y->divTest.BigIntModK N (x+y))
    printf "%A" (finalAns )
    0 // return an integer exit code
