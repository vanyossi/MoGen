/*
 *  Copyright (c) 2018 Iván Yossi <ghevan@gmail.com>
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 3 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <MacTypes.h>
#include "mogen_mop.h"

int main(int argc, char const *argv[]) {

    mop_print(mogen_mop(0, MOP_BIN | MOP_RESTRICTED));
    mop_print(mogen_mop(0, MOP_REAL | MOP_CONTIGUOUS));
    mop_print(mogen_mop(0, MOP_DYNAMIC | MOP_BIN | MOP_REAL));
    mop_print(mogen_mop(0, MOP_RESTRICTED | MOP_REAL | MOP_DYNAMIC));

    return 0;
}
