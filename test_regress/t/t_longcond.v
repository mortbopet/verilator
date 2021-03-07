// DESCRIPTION: Verilator: Verilog Test module
//
// This file ONLY is placed under the Creative Commons Public Domain, for
// any use, without warranty, 2019 by Wilson Snyder.
// SPDX-License-Identifier: CC0-1.0

module t ( /*AUTOARG*/
    // Inputs
    input wire clk,
    input wire inCond,
    input wire [31:0] inComp,
    output reg [31:0] out
    );

    integer i;
    reg [9:0] o_reg;
    reg o;

    reg [31:0] value;
    always @(posedge clk ) begin

        for (i=0; i < 10; i = i + 1) begin
            o = o ^ o_reg[i] ^ inCond;
        end

        if(o) begin
            value += inComp;
            value *= value;
        end else begin
            value -= inComp;
            value *= value;
        end
        out <= value;
        o_reg = o_reg << 1;
        o_reg[0] = o;
    end
endmodule
