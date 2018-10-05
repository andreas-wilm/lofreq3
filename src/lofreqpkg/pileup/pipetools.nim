## Provides several useful procedures for procedure chainging (making
## pipelines).
## IMPORTANT CAVEAT: Nim does not provide a way to ensure functional purity.
## Thus, all involved procedures could, in mathematical sense, have side
## effects. 

func compose*[A, B, C](f: func (x: A) : B,
                       g: func (y: B): C): (func (z: A): C) =
  ## Composes two procedures. In mathematical notation, it takes two procedures
  ## with the following type signatures:
  ## f: A -> B
  ## g: B -> C
  ## and returns a new procedure with a type signature:
  ## (f ∘ g): A -> C
  ## Where the resulting value (and side effects) of (f ∘ g)(x) are equal to
  ## those of f(g(x)).
  ## IMPORTANT CAVEAT: Nim does not provide a way to ensure functional purity.
  ## Thus, all given procedure could, in mathematical sense, have side effects. 
  return func (value: A): C = f(g(value))


func compose*[A, B](f: func (x: A) : B,
                    g: proc (y: B): void): (proc (z: A): void) =
  ## Composes two procedures. In mathematical notation, it takes two procedures
  ## with the following type signatures:
  ## f: A -> B
  ## g: B -> void 
  ## and returns a new procedure with a type signature:
  ## (f ∘ g): A -> void 
  ## Where the effects of (f ∘ g)(x) are equal to those of (f(g(x)). In both
  ## cases, nothing is returned.
  ## IMPORTANT CAVEAT: Nim does not provide a way to ensure functional purity.
  ## Thus, all given procedure could, in mathematical sense, have side effects.
  ## And the second one, given that it has a return value of void, most likely
  ## will. This overload exists because 'void' is not a type and the regular
  ## 'compose' procedure cannot be used.
  return proc (value: A): void = f(g(value))


func then*[A, B, C](f: func (x: A) : B, g: func (y: B): C): (func (z: A): C) =
  ## Composes two procedures in reverse order. It is equivalent to calling the
  ## compose funtion with flipped arguments. 
  return func (value: A): C = g(f(value))


func then*[A, B](f: func (x: A) : B,
                 g: func (y: B): void): (func (z: A): void) =
  ## Composes two procedures in reverse order. It is equivalent to calling the
  ## compose funtion with flipped arguments. This overload exists because 'void'
  ## is not a type and the regular 'then' procedure cannot be used.
  return proc (value: A): void = g(f(value))


func thenDo*[A, B](f: func (x: A) : B,
                   g: func (y: B): void): (func (z: A): B) =
  ## Returns a procedure returning the value returned by the first procedure
  ## after applying the second procedure to the said value.
  ## IMPORTANT CAVEAT: Nim does not provide a way to ensure functional purity.
  ## Thus, all given procedure could, in mathematical sense, have side effects.
  ## And the second one, given that it has a return value of void, most likely
  ## will.
  return proc (value: A): B =
    let firstResult = f(value)
    g(firstResult)
    return firstResult


func done*[A, B](f: func (x: A) : B): (proc (z: A): void) =
  ## Wraps the discard operator to enable easier chaining. Returns a procedure
  ## providing the same functionality as the procedure given as the argument,
  ## but which discards the result and doesn't return anything.
  return proc (value: A): void = discard f(value)
 

