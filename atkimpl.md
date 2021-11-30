# SEAL
攻击成功 库代码修改：加速防御　防御：3.2 stddev
# HElib
攻击仅在无同态操作时条件成功 库代码修改：暴露接口/绕过防御 防御：自带
# HEAAN
攻击成功 库代码修改：无 防御：3.2 stddev
# FullRNS-HEAAN
攻击条件成功 库代码修改：decode逻辑 防御：3.2 stddev
the original version contains the conversion of the first prime into 'double' type
this will introduce an error of magnitude exactly
(consider the bit pattern of the first prime:
1x...x0...01, with 62 bits in total and k-1 successive 0-bits ahead of the lowest 1-bit
the highest 53 bits are kept, and the lower 9 bits are discarded
which have the pattern 0...01(since k-1 is typically greater than 8)
BUT, the coeff is also converted to double, then subtracted by the first prime(this conversion itself is accurate)
following the rules for floating point arithmetic, the operand with a smaller exponent is scaled up so that its
exponent matches the other,
during this process, the lower 9 bits of the integral part are discarded,
and a literally huge error is introduced, since those bits can take any value
# PALISADE
攻击条件成功 库代码修改：绕过防御 防御：自带
# lattigo
攻击成功 库代码修改：暴露接口 防御：3.2 stddev