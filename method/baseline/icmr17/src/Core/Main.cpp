/* <Main.cpp>
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 ***********************************************************************************
 *
 * This file is part of Unsupervised Distance Learning Framework (UDLF).
 *
 * UDLF is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * UDLF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with UDLF.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <string>

#include "Exec.hpp"

int main(int argc, char* argv[]) {

    std::string filename;
    std::string method;
    std::string dataset_size;
    std::cout<<"process begin with"<< argc<<std::endl;
    if (argc == 4){

        filename = argv[1];
        method = argv[2];
	dataset_size  = argv[3];
    }else{
        std::cout<<"not enough input para"<<std::endl;
        return(4);
    }

    Exec& exec = Exec::getInstance();
    if (!exec.parseFile(filename, method, dataset_size)) {
        return 1;
    }
    exec.run();

    return 0;
}
