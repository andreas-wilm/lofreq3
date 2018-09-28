func compose*[A, B, C](f: func (x: A) : B, g: func (y: B): C): (func (z: A): C) =
  return func (value: A): C = f(g(value))


func composeImpure*[A, B](f: func (x: A) : B, g: proc (y: B): void): (proc (z: A): void) =
  return proc (value: A): void = f(g(value))


func then*[A, B, C](f: func (x: A) : B, g: func (y: B): C): (func (z: A): C) =
  return func (value: A): C = g(f(value))


func then*[A, B](f: func (x: A) : B, g: func (y: B): void): (func (z: A): void) =
  return proc (value: A): void = g(f(value))


func thenDo*[A, B](f: func (x: A) : B, g: func (y: B): void): (func (z: A): B) =
  return proc (value: A): B =
    let firstResult = f(value)
    g(firstResult)
    return firstResult


func done*[A, B](f: func (x: A) : B):(proc (z: A): void) =
  return proc (value: A): void = discard f(value)
 

#func compose*[A, B, C](f: func (x: A) : B, g: func (y: B): C): (func (z: A): C) =
#  return func (value: A): C = f(g(value))
#
#func composeImpure*[A, B](f: func (x: A) : B, g: proc (y: B): void): (proc (z: A): void) =
#  return proc (value: A): void = f(g(value))
#
#func then*[A, B, C](f: func (x: A) : B, g: func (y: B): C): (func (z: A): C) =
#  return func (value: A): C = g(f(value))
#
#func then*[A, B](f: func (x: A) : B, g: func (y: B): void): (func (z: A): B) =
#  return proc (value: A): B =
#    let firstResult = f(value)
#    g(firstResult)
#    return firstResult
#
#func finallyDo*[A, B](f: func (x: A) : B, g: proc (y: B): void): (proc (z: A): void) =
#  return proc (value: A): void = g(f(value))
#
#func finallyDo*[A, B, C](f: func (x: A) : B, g: func (y: B): C): (proc (z: A): void) =
#  return proc (value: A): void = discard g(f(value))
 
