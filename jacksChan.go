package main

import (
	"fmt"
	"github.com/wiless/vlib"
	"math"
)

var matlab *vlib.Matlab

type JakesModel struct {
	fs, fd       float64
	ts           float64
	beta, phi, f vlib.VectorF
	M, N         float64
}

type MultiTapFading []JakesModel

func NewMultiTapFading(N int, fs, fd float64) MultiTapFading {
	m := make([]JakesModel, N)

	for i := 0; i < N; i++ {
		m[i].Init(fs, fd)
	}
	return m
}

func (m *MultiTapFading) Generate(indx float64) vlib.VectorC {
	result := vlib.NewVectorC(len(*m))

	for i := 0; i < len(*m); i++ {
		result[i] = (*m)[i].Generate(indx)
	}

	return result
}

func main() {
	matlab = vlib.NewMatlab("test.m")
	channelVal := vlib.NewVectorC(1000)
	hn := NewMultiTapFading(5, 2000, 15)

	for t := 0; t < 10000; t++ {
		fmt.Println(t, " <=>  ", hn.Generate(float64(t)))
	}
	// for i := 0; i < channelVal.Size(); i++ {
	// 	channelVal[i] = model.Generate(float64(i))
	// }
	// rms := vlib.Norm2C(channelVal)
	// matlab.Export("h", channelVal.Scale(1/rms))
	// matlab.Command("close all;plot(abs(h))")
	matlab.Close()
	fmt.Println("..bye bye")
	fmt.Println(channelVal)
}

func (j *JakesModel) Init(fs, fd float64) {
	j.M = 20.0
	j.N = 4*j.M + 20
	j.fs = fs
	j.ts = 1 / fs
	j.fd = fd

	j.beta = vlib.NewVectorF(int(j.M))
	j.f = vlib.NewVectorF(int(j.M))
	for i := 0; i < j.beta.Size(); i++ {
		count := math.Pi * float64(i+1)
		j.beta[i] = count / (j.M + 1)
		j.f[i] = j.fd * math.Cos(2*count/j.N) // Set values of f(k,n)
	}

	x := vlib.VectorF(vlib.RandUFVec(len(j.beta) + 1))
	j.phi = x.Scale(2 * math.Pi).Add(-math.Pi)

}

func (j *JakesModel) Generate(m float64) complex128 {
	var zR, zI float64
	pi := math.Pi
	for i := 0; i < j.beta.Size(); i++ {

		zR += 2 * math.Cos(j.beta[i]) * math.Cos((2*pi*j.f[i]*m*j.ts)+j.phi[i])
		zI += 2 * math.Sin(j.beta[i]) * math.Cos((2*pi*j.f[i]*m*j.ts)+j.phi[i])
	}
	z1 := complex(zR, zI)
	z2 := complex(math.Sqrt(2)*math.Cos((2*math.Pi*j.fd*m*j.ts)+j.phi[j.beta.Size()]), math.Sqrt(2)*math.Sin((2*math.Pi*j.fd*m*j.ts)+j.phi[j.beta.Size()]))

	z := (z1 + z2) * complex(2.5/math.Sqrt(j.N), 0)
	return z
}

// func jacksm(samplingRate float64, m float64, fd float64, beta vlib.NewVectorF, phi vlib.NewVectorF) (complex128, vlib.NewVectorF, vlib.NewVectorF) {

// 	// samplingRate := 2000.0
// 	// m := 6.0

// 	// fmt.Println(sampling_period)
// 	// phi(n) ~ U[-pi,pi) is a RV
// 	//for m = 0:no_of_samples-1
// 	// for i:=1
// 	var z_term_sum_real float64
// 	var z_term_sum_imag float64
// 	for i := 0; i < beta.Size(); i++ {
// 		z_term_sum_real += 2 * math.Cos(beta[i]) * math.Cos((2*pi*f[i]*m*sampling_period)+phi[i])
// 		z_term_sum_imag += 2 * math.Sin(beta[i]) * math.Cos((2*pi*f[i]*m*sampling_period)+phi[i])
// 	}

// 	z_term_2_real := math.Sqrt(2) * math.Cos((2*pi*fd*m*sampling_period)+phi[beta.Size()])
// 	z_term_2_imag := math.Sqrt(2) * math.Sin((2*pi*fd*m*sampling_period)+phi[beta.Size()])

// 	z_real := (2 / math.Sqrt(N)) * (z_term_sum_real + z_term_2_real)
// 	z_imag := (2 / math.Sqrt(N)) * (z_term_sum_imag + z_term_2_imag)

// 	z := complex(z_real, z_imag)
// 	return z, beta, phi
// 	// fmt.Println(z)
// 	//end
// 	// rms_value = sqrt(((norm(abs(z)))^2)/no_of_samples);
// 	// z = z/rms_value ;       % Normalisation is being done
// }
