/**
MINMAXGD
Copyright (C) {{ 2007 }}  {{ Univeristy of Angers, L. Hardouin }}

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "..\include\matrixCOO.h"

namespace mmgd
{
	void matrixCOO::affecte(matrixCOO &list)
	{
		int i;

		this->row = list.row;
		this->col = list.col;
		this->size = list.size;
		if (this->capacity < list.size)
		{
			capacity = list.capacity;

			if (this->Data != NULL)
			{
				delete[] Data;
			}

			Data = new Ttriplet[capacity];
		}
		if (this->Data == NULL)
		{
			mem_limite l(15);
			throw(l);
		}

		for (i = 0; i < size; i++)
		{
			this->Data[i].s = list.Data[i].s;

			Data[i].i = list.Data[i].i;
			Data[i].j = list.Data[i].j;
		}
	}

	// constructeur par defaut de la classe smatrix (matrice 1x1 contenant
	// epsilon)
	matrixCOO::matrixCOO(void)
	{
		row = 0;
		col = 0;
		size = 0;
		capacity = 0;

		this->Data = NULL;
	}

	matrixCOO::matrixCOO(int i, int j) // constructeur initialisant
	{

		row = i;
		col = j;
		size = 0;
		capacity = 0; // par défaut on prend une matrice de capcité 25%
		this->Data = NULL;
	}

	matrixCOO::matrixCOO(const matrixCOO &a) // constructeur initialisation par une autre matrice
	{

		affecte((matrixCOO &)a);
	}

	matrixCOO::matrixCOO(const serie &a) // constructeur initialisation par une serie
	{
	}

	matrixCOO::matrixCOO(poly &a) // constructeur initialisation par un polyn�e
	{
	}

	matrixCOO::matrixCOO(gd &a) // constructeur initialisation par un monome
	{
	}

	matrixCOO::~matrixCOO() // destructeur
	{
		int n;

		delete[] Data;

		row = 0;
		col = 0;
		size = 0;
		capacity = 0;
		Data = NULL;
	}

	int matrixCOO::matrixCOO::geti(int k)
	{
		return Data[k].i;
	}
	int matrixCOO::getj(int k)
	{
		return Data[k].j;
	}
	int matrixCOO::getcapacity()
	{
		return capacity;
	}
	int matrixCOO::comp(int ai, int aj, int bi, int bj)
	{

		if (ai > bi)
			return 1;
		else
		{
			if (ai < bi)
				return -1;
			else
			{
				if (aj > bj)
					return 1;
				else
				{
					if (aj < bj)
						return -1;
					else
						return 0;
				}
			}
		}
	}
	void matrixCOO::pushback(Ttriplet &newdata)
	{
		Ttriplet *newvect;
		int k;
		int oldsize = size;
		size++;
		if (size > capacity)
		{
			capacity = capacity + (row * col);
			// cout<<"capacity"<<capacity<<endl;
			newvect = new Ttriplet[capacity];

			if (newvect == NULL)
			{
				mem_limite l(15);
				throw(l);
			}

			for (k = 0; k < oldsize; k++)
				newvect[k] = Data[k];
			delete[] Data;
			Data = newvect;
		}

		Data[size - 1].j = newdata.j;
		Data[size - 1].i = newdata.i;
		Data[size - 1].s = newdata.s;
	}
	int matrixCOO::find(int i, int j)
	{
		int imax = size - 1;
		int imin = 0;
		int imid;
		while (imax >= imin)
		{
			/* calculate the midpoint for roughly equal partition */
			imid = (imin + imax) / 2;

			// determine which subarray to search
			if (Data[imid].i < i)
				// change min index to search upper subarray
				imin = imid + 1;
			else if (Data[imid].i > i)
				// change max index to search lower subarray
				imax = imid - 1;
			else
			// we have to look for column
			{
				if (Data[imid].j < j)
					// change min index to search upper subarray
					imin = imid + 1;
				else if (Data[imid].j > j)
					imax = imid - 1;
				else
				// we have found
				{
					return imid;
				}
			}
		}
		// key not found
		return -1;
	}
	matrixCOO &matrixCOO::operator=(const matrixCOO &a)
	// initialise avec un objet smatrixCOO, surdefinition du =
	{

		if (&a == this)
			return *this; // si a est = de la matrice courante

		affecte((matrixCOO &)a);
		return *this;
	}
	matrixCOO &matrixCOO::operator=(smatrix &a)
	// initialise avec une matrice de serie, la matrice destination sera détruite si elle existait déja
	{
		serie eps;
		Ttriplet aij;
		int i, j, k;

		row = a.getrow();
		col = a.getcol();
		size = 0;

		for (i = 0; i < a.getrow(); i++)
		{
			for (j = 0; j < getcol(); j++)
			{
				if (!(a(i, j) == eps))
				{
					aij.s = a(i, j);
					aij.i = i;
					aij.j = j;
					pushback(aij);
				}
			}
		}
	}

	matrixCOO &matrixCOO::operator=(serie &a) // surdefinition du =, permet d'initialiser avec une serie cast serei matrice
	{
	}

	matrixCOO &matrixCOO::operator=(poly &p1) // initialise avec un polynome cast polynome->matrice
	{
	}

	matrixCOO &matrixCOO::operator=(gd &gd1) // initialise avec un monome cast monome->matrice
	{
	}

	serie &matrixCOO::operator()(int i, int j)
	{ // serie eps;
		int index;
		if (i >= row || j >= col || i < 0 || j < 0)
		{
			taille_incorrecte number(1);
			throw(number);
		}

		// Recherche Dichotomique
		index = find(i, j);
		if (index == -1)
		// it doesn't exist, we have to inseret a series, but the problem now, we have some epsilon in the matrix
		{
			return (eps);
		}

		return Data[index].s;
	}

	int matrixCOO::operator==(const matrixCOO &M)
	{
		int k;

		if (M.col != this->col || M.row != this->row)
		{
			taille_incorrecte number(1);
			throw(number);
		}

		if (M.size != this->size)
			return 0;

		for (k = 0; k != this->size; k++)
		{
			if (this->Data[k].i != M.Data[k].i || this->Data[k].j != M.Data[k].j)
			{

				return 0;
			}

			if (!(this->Data[k].s == M.Data[k].s))
			{
				return 0;
			}
		}

		return 1;
	}

	int matrixCOO::operator!=(const matrixCOO &M)
	{
		int k;

		if (M.col != this->col || M.row != this->row)
		{
			taille_incorrecte number(1);
			throw(number);
		}

		if (M.size != this->size)
			return 1;

		for (k = 0; k != this->size; k++)
		{
			if (this->Data[k].i != M.Data[k].i || this->Data[k].j != M.Data[k].j)
			{
				return 1;
			}
			if (!(this->Data[k].s == M.Data[k].s))
			{
				return 1;
			}
		}

		return 0;
	}

	std::ostream &operator<<(std::ostream &flot, matrixCOO &a)
	// surdefinition de <<
	// affichage
	{
		int k;

		for (k = 0; k < a.size; k++)
		{

			flot << "[" << a.Data[k].i << "," << a.Data[k].j << "] = " << a.Data[k].s << "\n";
		}

		return flot;
	}

	std::fstream &operator<<(std::fstream &flot, matrixCOO &a)
	// surdefinition de <<
	// pour fichier texte
	{
		int k;

		for (k = 0; k < a.size; k++)
		{

			flot << "[" << a.Data[k].i << "," << a.Data[k].j << "] = " << a.Data[k].s << "\n";
		}

		return flot;
	}

	matrixCOO oplus(matrixCOO &a, matrixCOO &b)
	{
		int i, j, index;
		Ttriplet temp;
		matrixCOO result;
		result.row = a.getrow();
		result.col = a.getcol();
		result.size = 0;

		i = 0;
		j = 0;
		if ((a.col == b.col) && (a.row == b.row))
		{

			while ((i != a.size) && (j != b.size))
			{
				index = a.comp(a.geti(i), a.getj(i), b.geti(j), b.getj(j));

				if (index == 0)
				{
					temp.s = oplus(a.Data[i].s, b.Data[j].s);
					temp.i = a.geti(i);
					temp.j = a.getj(i);

					result.pushback(temp);

					i++;
					j++;
				}
				else
				{
					if (index == -1)
					{
						temp.s = a.Data[i].s;
						temp.i = a.geti(i);
						temp.j = a.getj(i);
						result.pushback(temp);

						i++;
					}
					else
					{

						temp.s = b.Data[j].s;
						temp.i = b.geti(j);
						temp.j = b.getj(j);
						result.pushback(temp);

						j++;
					}
				}
			}
			while (i != a.size) // il ne reste plus que des A
			{
				temp.s = a.Data[i].s;
				temp.i = a.geti(i);
				temp.j = a.getj(i);
				result.pushback(temp);

				i++;
			}
			while (j != b.size) // il ne reste plus que des B
			{

				temp.s = b.Data[j].s;
				temp.i = b.geti(j);
				temp.j = b.getj(j);
				result.pushback(temp);

				j++;
			}

			return (result);
		}
		else
		{
			mem_limite l(27);
			throw(l);
		}
	}

	matrixCOO inf(matrixCOO &a, matrixCOO &b)
	{
		int i = 0, j = 0, index = 0;
		Ttriplet temp;
		matrixCOO result;
		result.row = a.getrow();
		result.col = a.getcol();
		result.size = 0;
		if ((a.col == b.col) && (a.row == b.row))
		{
			while ((i != a.size) && (j != b.size))
			{
				index = a.comp(a.geti(i), a.getj(i), b.geti(j), b.getj(j));

				if (index == 0)
				{
					temp.s = inf(a.Data[i].s, b.Data[j].s);
					temp.i = a.geti(i);
					temp.j = a.getj(i);

					result.pushback(temp);

					i++;
					j++;
				}
				else
				{
					if (index == -1)
					{

						i++;
					}
					else
					{

						j++;
					}
				}
			}

			return (result);
		}
		else
		{
			mem_limite l(27);
			throw(l);
		}
	}
	void transpose(matrixCOO &a, matrixCOO &b)
	{
		int n = 0, i = 0, j = 0, currentb = 0;
		n = a.getsize();
		b.row = a.getcol();
		b.col = a.getrow();
		serie seps;
		seps.init(epsilon, epsilon, e);
		Ttriplet temp;
		b.size = n;
		if (n > 0)
		{
			currentb = 1;
			for (i = 0; i < a.col; i++)
			{
				for (j = 1; j <= n; j++)
				{
					if (a.Data[j].j == i)
					{
						temp.i = a.getj(j);
						temp.j = a.geti(j);
						temp.s = a.Data[j].s;
						if (!((temp.s == seps) || (a(i, j) == seps)))
						{

							b.pushback(temp);
						}
						currentb++;
					}
				}
			}
		}
	}

	void fast_transpose(matrixCOO &a, matrixCOO &b)
	{

		int row_terms[MAX_COL], starting_pos[MAX_COL];
		int i = 0, j = 0, num_cols = a.col, num_terms = a.size;
		b.row = num_cols;
		b.col = a.row;
		b.size = num_terms;
		if (num_terms > 0)
		{
			for (i = 0; i < num_cols; i++)
			{
				row_terms[i] = 0;
			}

			for (i = 0; i < num_terms; i++)
			{
				row_terms[a.Data[i].j]++;
				starting_pos[0] = 1;
			}

			for (i = 0; i < num_cols; i++)
			{
				starting_pos[i] = starting_pos[i - 1] + row_terms[i - 1];
			}
			for (i = 0; i < num_terms; i++)
			{
				j = starting_pos[a.Data[i].j]++;
				b.Data[j].i = a.Data[i].j;
				b.Data[j].j = a.Data[i].i;
				b.Data[j].s = a.Data[i].s;
			}
		}
	}

} // fin namespace
