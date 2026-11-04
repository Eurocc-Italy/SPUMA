/*---------------------------------------------------------------------------*\

License
    This file is part of the TurboWG project:
    https://sourceforge.net/projects/turbowg

    The files from the TurboWG project are free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    This source code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this source code.  If not, see <http://www.gnu.org/licenses/>.

    Special notice:
    This file is part of a port of the mixingPlane interface source code from
    the foam-extend project (https://sourceforge.net/projects/foam-extend/) to
    OpenFOAM (https://www.openfoam.com/). Most of the original source code has
    been preserved, including the comments and authorship/IP notices.

Contributors/Copyright:
    Martin Beaudoin, (2023)

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "twg_CSV.H"
#include "DynamicList.H"
#include "ListOps.H"

#define NO_REF_COLUMN -1

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::labelList Foam::Function1Types::twg_CSV<Type>::getComponentColumns
(
    const word& name,
    const dictionary& dict
)
{
    // Writing of columns was forced to be ASCII,
    // do the same when reading

    labelList cols;
    ITstream& is = dict.lookup(name);
    is.format(IOstreamOption::ASCII);
    is >> cols;
    dict.checkITstream(is, name);

    if (cols.size() != pTraits<Type>::nComponents)
    {
        FatalIOErrorInFunction(dict)
            << name << " with " << cols
            << " does not have the expected length "
            << pTraits<Type>::nComponents << nl
            << exit(FatalIOError);
    }
    return cols;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::twg_CSV<Type>::readValue
(
    const UList<string>& strings
) const
{
    Type result;
    if (std::is_integral<Type>::value)
    {
        // nComponents == 1
        setComponent(result, 0) = readLabel(strings[componentColumns_[0]]);
    }
    else
    {
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
        {
            setComponent(result, cmpt) =
                readScalar(strings[componentColumns_[cmpt]]);
        }
    }
    return result;
}


template<class Type>
void Foam::Function1Types::twg_CSV<Type>::read()
{
    fileName expandedFile(fName_);
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(expandedFile.expand()));
    ISstream& is = isPtr();

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open twg_CSV file for reading."
            << exit(FatalIOError);
    }

    const label maxEntry =
        max(refColumn_, componentColumns_[findMax(componentColumns_)]);

    string line;
    label lineNo = 0;

    // Skip header
    for (label i = 0; i < nHeaderLine_; ++i)
    {
        is.getLine(nullptr);
        ++lineNo;
    }

    DynamicList<Tuple2<scalar, Type>> values;
    DynamicList<string> strings(maxEntry+1);  // reserve

    while (is.good())
    {
        is.getLine(line);
        ++lineNo;

        // Skip commented out lines
        if
        (
            !line.empty() && line.at(0) == comment_
        )
        {
            continue;
        }

        strings.clear();

        std::size_t pos = 0;

        for
        (
            label n = 0;
            (pos != std::string::npos) && (n <= maxEntry);
            ++n
        )
        {
            if (mergeSeparators_)
            {
                bool found = false;
                while (!found)
                {
                    const auto nPos = line.find(separator_, pos);

                    if ((nPos != std::string::npos) && (nPos - pos == 0))
                    {
                        pos = nPos + 1;
                    }
                    else
                    {
                        found = true;
                    }
                }
            }

            const auto nPos = line.find(separator_, pos);

            if (nPos == std::string::npos)
            {
                strings.push_back(line.substr(pos));
                pos = nPos;
            }
            else
            {
                strings.push_back(line.substr(pos, nPos - pos));
                pos = nPos + 1;
            }
        }

        if (strings.size() <= 1)
        {
            break;
        }

        if (strings.size() <= maxEntry)
        {
            FatalErrorInFunction
                << "Not enough columns near line " << lineNo
                << ". Require " << (maxEntry+1) << " but found "
                << strings << nl
                << exit(FatalError);
        }

        // Use refColumn_ values as increasing time values, or if
        // refColumn_ == NO_REF_COLUMN, replace by values based on
        // the current line number of the CSV file
        scalar x = (refColumn_ != NO_REF_COLUMN) ?
            readScalar(strings[refColumn_]) :           // As specified by user
            scalar(is.lineNumber() - nHeaderLine_ - 2); // Ramp starting at 0

        Type value = readValue(strings);

#if OPENFOAM < 2302
        values.append(Tuple2<scalar,Type>(x, value));
#else
        values.emplace_back(x, value);
#endif
    }

    this->table_.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::twg_CSV<Type>::twg_CSV
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr,
    const fileName& fName
)
:
    TableBase<Type>(entryName, dict, obrPtr),
    nHeaderLine_(dict.getOrDefault<label>("nHeaderLine", 0)),
    refColumn_(dict.getOrDefault<label>("refColumn", NO_REF_COLUMN)),
    componentColumns_(getComponentColumns("componentColumns", dict)),
    separator_(dict.getOrDefault<string>("separator", ",")[0]),
    comment_(dict.getOrDefault<string>("commentChar", "#")[0]),
    mergeSeparators_(dict.getOrDefault<bool>("mergeSeparators", false)),
    fName_(fName.empty() ? dict.get<fileName>("file") : fName)
{
    // Some sanity checks
    if(separator_ == comment_)
    {
        WarningInFunction
            << "The CSV field separator is the same as the comment character."
            << " Disabling the parsing of comments."
            << nl << endl;

        comment_ = '\0';
    }

    read();

    TableBase<Type>::initialise();
}


template<class Type>
Foam::Function1Types::twg_CSV<Type>::twg_CSV(const twg_CSV<Type>& csv)
:
    TableBase<Type>(csv),
    nHeaderLine_(csv.nHeaderLine_),
    refColumn_(csv.refColumn_),
    componentColumns_(csv.componentColumns_),
    separator_(csv.separator_),
    comment_(csv.comment_),
    mergeSeparators_(csv.mergeSeparators_),
    fName_(csv.fName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fileName& Foam::Function1Types::twg_CSV<Type>::fName() const
{
    return fName_;
}


template<class Type>
void Foam::Function1Types::twg_CSV<Type>::writeEntries(Ostream& os) const
{
    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    os.writeEntry("nHeaderLine", nHeaderLine_);

    if(refColumn_ != NO_REF_COLUMN)
    {
        os.writeEntry("refColumn", refColumn_);
    }

    // Force writing labelList in ASCII
    const auto oldFmt = os.format(IOstreamOption::ASCII);
    os.writeEntry("componentColumns", componentColumns_);
    os.format(oldFmt);

    os.writeEntry("separator", string(separator_));
    os.writeEntry("commentChar", string(comment_));
    os.writeEntry("mergeSeparators", mergeSeparators_);
    os.writeEntry("file", fName_);
}


template<class Type>
void Foam::Function1Types::twg_CSV<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //
