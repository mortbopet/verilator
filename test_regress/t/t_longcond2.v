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
    output reg [31:0] out1,
    output reg [31:0] out2
    );

    integer i;
    reg [31:0] v1, v2;
    reg [9:0] o_reg;
    reg o, o_end;

    always @(*) begin 
        o_end = o_reg[9] | inCond;
    end


    always @(*) begin 
        for (i=0; i < 10; i = i + 1) begin
            o = o ^ o_reg[i] ^ inCond;
        end
    end

    always @(*) begin
        v1 = 0;
        if(o_end) begin
            for (i=0; i < 10; i = i + 1) begin
                v1 += inComp;
                v1 *= v1;
            end
        end
    end

    always @(*) begin
        v2 = 0;
        if(o_end) begin
            for (i=0; i < 10; i = i + 1) begin
                v2 += inComp;
                v2 *= v2;
            end
        end
    end

    always @(posedge clk ) begin
        out1 <= v1;
        out2 <= v2;
        o_reg = o_reg << 1;
        o_reg[0] = o;
    end
endmodule
