/*
 * lookup_table.h
 *
 *  Created on: Jun 29, 2015
 *      Author: IEUser
 */

#ifndef LOOKUP_TABLE_H_
#define LOOKUP_TABLE_H_

#define W 1000000
//#define W 2000000

//#define C1  W / (2 * PI)
#define C1 159154.943091895
//#define C2 (2 * PI) / W
#define C2 0.0000062831853967

//float sin_lut(float);
//float cos_lut(float);
//float sinc_lut(float);

double sin_lut(double);
double cos_lut(double);
double sinc_lut(double);


#endif /* LOOKUP_TABLE_H_ */
