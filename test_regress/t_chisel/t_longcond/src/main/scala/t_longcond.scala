import chisel3._
import chisel3.Driver
import chisel3.util.ShiftRegister
import chisel3.util.Cat

class T_longcond(depth : Int) extends Module {
  val nFuncUnits = (math.pow(2, depth) - 1).toInt
  val io = IO(new Bundle {
    val in1 = Input(UInt(32.W))
    val in2 = Input(UInt(32.W))
    val doComp = Input(UInt(nFuncUnits.W))
    val out = Output(UInt(32.W))
  })

  val COND_COMPLEXITY = 10 // defines the computation size of the comparison value
  val COMP_COMPLEXITY = 10 // defines the computation size of the computation value

  val regs = Reg(Vec(nFuncUnits, (UInt(32.W))))
  val adders = Wire(Vec(nFuncUnits, UInt(32.W)))
  val addCond = Wire(Vec(nFuncUnits, Bool()))
  val regEn = Wire(Bool())

  io.out := adders.slice((nFuncUnits + 1)/2, nFuncUnits).reduce(_*_)
  regs := adders.zipWithIndex.map({case (adder, i) => Mux(regEn, adder, regs(i))})
  regEn := addCond.reduce(_ & _)
  addCond := (heavyBooleanComp(io.in1, io.in2, regs(0))) +: (1 until nFuncUnits).map{i => heavyBooleanComp(regs(i), io.in2, regs(i-1))}

  for(i <- 0 until nFuncUnits) {
    val preIdx = (i-1)/2 
    val op1 = if (i == 0) io.in1 else adders(preIdx)
    adders(i) := Mux(addCond(i), op1 + regs(i), regs(i) * 2.U)
  }

  private def heavyBooleanComp(i: UInt, j: UInt, k: UInt) : Bool = {
    def compStuff(b: UInt, c: UInt, d: UInt, round: Int) : UInt = {
      round % 4 match {
        case 0 => (b & c) | (~b & d)
        case 1 => (b & d) | (c & ~d)
        case 2 => (b ^ c ^ d)
        case 3 => c ^ (b | ~d)
      }
    }

    (0 until COND_COMPLEXITY).foldLeft(i){case (acc, i) => compStuff(acc, j, k, i)} === 42.U
  }

  // private def heavyCombComp()
}


object T_longcond extends App {
  (new chisel3.stage.ChiselStage).emitVerilog(new T_longcond(2))
}