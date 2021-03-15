import chisel3._
import chisel3.Driver
import chisel3.util.ShiftRegister
import chisel3.util.Cat

class T_longcond(depth : Int) extends Module {
  val nFuncUnits = (math.pow(2, depth) - 1).toInt
  val io = IO(new Bundle {
    val in = Input(UInt(32.W))
    val doComp = Input(UInt(nFuncUnits.W))
    val out = Output(UInt(32.W))
  })

  val regs = Reg(Vec(nFuncUnits, (UInt(32.W))))
  val adders = Wire(Vec(nFuncUnits, UInt(32.W)))
  val addCond = Wire(Vec(nFuncUnits, Bool()))
  val regEn = Wire(Bool())

  io.out := adders.slice((nFuncUnits + 1)/2, nFuncUnits).reduce(_*_)
  regs := adders.zipWithIndex.map({case (adder, i) => Mux(regEn, adder, regs(i))})
  regEn := addCond.reduce(_ & _)
  addCond := (io.in > regs(0)) +: (1 until nFuncUnits).map{i => regs(i) > regs(i-1)}

  for(i <- 0 until nFuncUnits) {
    val preIdx = (i-1)/2 
    val op1 = if (i == 0) io.in else adders(preIdx)
    adders(i) := Mux(addCond(i), op1 + regs(i), regs(i) * 2.U)
  }
}


object T_longcond extends App {
  (new chisel3.stage.ChiselStage).emitVerilog(new T_longcond(3))
}