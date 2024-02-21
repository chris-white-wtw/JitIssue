# JitIssue

Demonstrates what appears to be a bug in JIT code generation in .NET 8.

## Usage

From the project/solution directory:

```
dotnet run -c release
```

## Observed output

When running with default settings, output similar to the following is seen:

```
Exception on iteration 93732, entry LogEntry { Ticks = 511308229, ThreadId = 298, aa = 1, bb = 2, tolerance = 0, yy0 = 0.5580284387721988 }
aplusb=3, a=1, b=2, a for power series=2
JitIssue.UserArgumentException: Both parameters for Beta function must be positive
   at JitIssue.CephesPort.IncompleteBetaIntegral(Double aa, Double bb, Tolerance tolerance, Double xx) in C:\git\JitIssue\CephesPort.cs:line 1668
   at JitIssue.CephesPort.InverseIncompleteBetaIntegral(Double aa, Double bb, Tolerance tolerance, Double yy0) in C:\git\JitIssue\CephesPort.cs:line 644
   at JitIssue.Program.Main(String[] args) in C:\git\JitIssue\Program.cs:line 15
```

When running with `DOTNET_TieredPGO=0`:

```
Exception on iteration 32967, entry LogEntry { Ticks = 510478075, ThreadId = 12, aa = 1, bb = 2, tolerance = 0, yy0 = 0.9860168987440667 }
aplusb=1, a=0, b=1, a for power series=2
JitIssue.UserArgumentException: Overflow(3) in Gamma function: x=0 on entry
   at JitIssue.CephesPort.Gamma(Double x) in C:\git\JitIssue\CephesPort.cs:line 1286
   at JitIssue.CephesPort.IncompleteBetaIntegralPowerSeries(Double a, Double b, Tolerance tolerance, Double x) in C:\git\JitIssue\CephesPort.cs:line 1581
   at JitIssue.CephesPort.IncompleteBetaIntegral(Double aa, Double bb, Tolerance tolerance, Double xx) in C:\git\JitIssue\CephesPort.cs:line 1697
   at JitIssue.CephesPort.InverseIncompleteBetaIntegral(Double aa, Double bb, Tolerance tolerance, Double yy0) in C:\git\JitIssue\CephesPort.cs:line 644
   at JitIssue.Program.Main(String[] args) in C:\git\JitIssue\Program.cs:line 15
```

Note:
* Running a debug build does not throw an exception.
* The issue does not occur if the target framework is changed to net6.0.

## Analysis

### Exception in IncompleteBetaIntegral

In the first case, the exception is thrown because of a failed check in `IncompleteBetaIntegral`.

At this point, the disassembly of the start of the method looks like this:

```
00007FF84D8D6CE0  push        rsi
00007FF84D8D6CE1  push        rbx  
00007FF84D8D6CE2  sub         rsp,148h  
00007FF84D8D6CE9  vzeroupper  
00007FF84D8D6CEC  vmovaps     xmmword ptr [rsp+130h],xmm6  
00007FF84D8D6CF5  vmovaps     xmmword ptr [rsp+120h],xmm7  
00007FF84D8D6CFE  vmovaps     xmmword ptr [rsp+110h],xmm8  
00007FF84D8D6D07  vmovaps     xmmword ptr [rsp+100h],xmm9  
00007FF84D8D6D10  vmovaps     xmmword ptr [rsp+0F0h],xmm10  
00007FF84D8D6D19  vmovaps     xmmword ptr [rsp+0E0h],xmm11  
00007FF84D8D6D22  vmovaps     xmmword ptr [rsp+0D0h],xmm12  
00007FF84D8D6D2B  vmovaps     xmmword ptr [rsp+0C0h],xmm13  
00007FF84D8D6D34  vmovaps     xmmword ptr [rsp+0B0h],xmm14  
00007FF84D8D6D3D  vmovq       xmm1,r8         ; bb := tolerance.Value  ??!
00007FF84D8D6D42  vmovaps     xmm6,xmm0       ; aa
00007FF84D8D6D46  vmovaps     xmm8,xmm1       ; tolerance.Value (because of ??! above)
00007FF84D8D6D4A  vmovaps     xmm7,xmm3       ; xx
00007FF84D8D6D4E  vxorps      xmm0,xmm0,xmm0  
00007FF84D8D6D52  vucomisd    xmm0,xmm6       ; 0 vs aa
00007FF84D8D6D56  jae         JitIssue.CephesPort.IncompleteBetaIntegral(Double, Double, Tolerance, Double)+05AAh (07FF84D8D728Ah)  
00007FF84D8D6D5C  vxorps      xmm0,xmm0,xmm0  
00007FF84D8D6D60  vucomisd    xmm0,xmm8       ; 0 vs tolerance.Value (should be bb!)
00007FF84D8D6D65  jae         JitIssue.CephesPort.IncompleteBetaIntegral(Double, Double, Tolerance, Double)+05AAh (07FF84D8D728Ah)  
```

Note the instruction `vmovq xmm1,r8`: this appears to be accidentally overwriting parameter `bb` with the value of
parameter `tolerance`.

For comparison, the disassembly when the method is first executed looks like this:

```
[bytes 55 48 which I imagine are push rsi, push rbx as above]
00007FF849B40B92  sub         esp,230h  
00007FF849B40B98  vzeroupper  
00007FF849B40B9B  lea         rbp,[rsp+230h]  
00007FF849B40BA3  xor         eax,eax  
00007FF849B40BA5  mov         qword ptr [rbp-208h],rax  
00007FF849B40BAC  vxorps      xmm4,xmm4,xmm4  
00007FF849B40BB0  vmovdqa     xmmword ptr [rbp-200h],xmm4  
00007FF849B40BB8  vmovdqa     xmmword ptr [rbp-1F0h],xmm4  
00007FF849B40BC0  mov         rax,0FFFFFFFFFFFFFE20h                                                                                                       ; loop for zeroing some stack
00007FF849B40BCA  vmovdqa     xmmword ptr [rax+rbp],xmm4  
00007FF849B40BCF  vmovdqa     xmmword ptr [rbp+rax+10h],xmm4  
00007FF849B40BD5  vmovdqa     xmmword ptr [rbp+rax+20h],xmm4  
00007FF849B40BDB  add         rax,30h  
00007FF849B40BDF  jne         JitIssue.CephesPort.IncompleteBetaIntegral(Double, Double, Tolerance, Double)+03Ah (07FF849B40BCAh)     ; end of stack-zeroing loop
00007FF849B40BE1  mov         qword ptr [rbp+20h],r8        ; tolerance.Value
00007FF849B40BE5  vmovsd      qword ptr [rbp+10h],xmm0      ; aa
00007FF849B40BEA  vmovsd      qword ptr [rbp+18h],xmm1      ; bb
00007FF849B40BEF  vmovsd      qword ptr [rbp+28h],xmm3      ; xx
00007FF849B40BF4  cmp         dword ptr [7FF849287188h],0  
00007FF849B40BFB  je          JitIssue.CephesPort.IncompleteBetaIntegral(Double, Double, Tolerance, Double)+072h (07FF849B40C02h)  
00007FF849B40BFD  call        00007FF8A79BDA80  
00007FF849B40C02  vxorps      xmm0,xmm0,xmm0  
00007FF849B40C06  vucomisd    xmm0,mmword ptr [rbp+10h]     ; 0 vs aa
00007FF849B40C0B  jae         JitIssue.CephesPort.IncompleteBetaIntegral(Double, Double, Tolerance, Double)+088h (07FF849B40C18h)  
00007FF849B40C0D  vxorps      xmm0,xmm0,xmm0  
00007FF849B40C11  vucomisd    xmm0,mmword ptr [rbp+18h]     ; 0 vs bb
00007FF849B40C16  jb          JitIssue.CephesPort.IncompleteBetaIntegral(Double, Double, Tolerance, Double)+0D9h (07FF849B40C69h)  
```

### Exception in Gamma

This appears to have a very similar cause to the first case.

`IncompleteBetaIntegralPowerSeries` is called with `a`=2, but when this in turn calls `Gamma(a)`, `a` is 0.

A disassembly of the start of `IncompleteBetaIntegralPowerSeries` again shows a `vmovq xmm0, r8` instruction overwriting
one parameter with another:

```
00007fff`9e228672 4881ecc8000000     sub     rsp, 0C8h
00007fff`9e228679 c5f877             vzeroupper 
00007fff`9e22867c c5f829b424b0000000 vmovaps xmmword ptr [rsp+0B0h], xmm6
00007fff`9e228685 c5f829bc24a0000000 vmovaps xmmword ptr [rsp+0A0h], xmm7
00007fff`9e22868e c57829842490000000 vmovaps xmmword ptr [rsp+90h], xmm8
00007fff`9e228697 c578298c2480000000 vmovaps xmmword ptr [rsp+80h], xmm9
00007fff`9e2286a0 488dac24d0000000   lea     rbp, [rsp+0D0h]
00007fff`9e2286a8 c4c1f96ec0         vmovq   xmm0, r8           ; looks wrong
00007fff`9e2286ad c5f828f0           vmovaps xmm6, xmm0
00007fff`9e2286b1 c5f828f9           vmovaps xmm7, xmm1
00007fff`9e2286b5 c57828c3           vmovaps xmm8, xmm3
```
