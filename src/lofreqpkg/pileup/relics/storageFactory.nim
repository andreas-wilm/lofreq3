import positionData
import iStorage
import slidingTable
import slidingDeque

proc getStorage*(implementation: string, size: int = 2000,
                 submitProc: proc (data: PositionData): void
                ): IStorage =
  case implementation
    of "deque": newSlidingDeque(size, submitProc).getIStorage()
    of "table": newSlidingTable(size, submitProc).getIStorage()
    else: raise newException(ValueError, 
      "Implementation '" & implementation & "' is not supported.")