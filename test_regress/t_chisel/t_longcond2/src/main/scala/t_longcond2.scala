import chisel3._
import chisel3.Driver
import chisel3.util.ShiftRegister
import chisel3.util.Cat

class T_longcond2(depth : Int) extends Module {
  val nFuncUnits = (math.pow(2, depth) - 1).toInt
  val io = IO(new Bundle {
    val in1 = Input(UInt(32.W))
    val out1 = Output(UInt(32.W))
    val out2 = Output(UInt(32.W))
  })

  val COND_COMPLEXITY = 3 // defines the computation size of the comparison value
  val COMP_COMPLEXITY = 3 // defines the computation size of the computation value

  
  val logReg = Reg(Bool())
  val log = Wire(Bool())
  log := io.in1.asBools().foldLeft(logReg){_ ^ _}
  logReg := log
  io.out1 := io.in1.asBools().foldLeft(0.S){case (acc, b) => 
    val vv = Wire(SInt())
    when (log ^ b) {vv :=1.S} .otherwise ({vv := -1.S })
    acc +io.in1*(vv)
  }.asUInt()
  io.out2 := io.in1.asBools().foldLeft(0.S){case (acc, b) =>
    val vv = Wire(SInt())
    when (log ^ b) {vv :=1.S} .otherwise ({vv := -1.S })
    acc +io.in1*(vv)
  }.asUInt()
}


object T_longcond2 extends App {
  (new chisel3.stage.ChiselStage).emitVerilog(new T_longcond2(3))
}