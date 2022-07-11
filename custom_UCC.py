# This code is part of Qiskit.
#
# (C) Copyright IBM 2021, 2022.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""
A custom Unitary Coupled-Cluster Ansatz.
"""

import logging
from functools import partial
from typing import Callable, List, Optional, Sequence, Tuple, Union

from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import EvolvedOperatorAnsatz
from qiskit.opflow import EvolutionBase, PauliTrotterEvolution

from qiskit_nature import QiskitNatureError
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.operators.second_quantization import FermionicOp, SecondQuantizedOp

from fermionic_excitation_generator import generate_fermionic_excitations


logger = logging.getLogger(__name__)


class custom_UCC(EvolvedOperatorAnsatz):
    r"""The Unitary Coupled-Cluster Ansatz. For more information, see [1].

    This Ansatz is an `EvolvedOperatorAnsatz` given by :math:`e^{T - T^{\dagger}}` where
    :math:`T` is the *cluster operator*. This cluster operator generally consists of excitation
    operators which are generated by
    :meth:`~qiskit_nature.circuit.library.ansatzes.utils.generate_fermionic_excitations`.
    This method constructs the requested excitations based on a `HartreeFock` reference state by
    default.
    You can also use a custom excitation generator method by passing a callable to `excitations`.

    A utility class :class:`UCCSD` exists, which is equivalent to:

    .. code-block:: python

        uccsd = UCC(excitations='sd', alpha_spin=True, beta_spin=True, max_spin_excitation=None)

    If you want to use a tailored Ansatz, you have multiple options to do so. Below, we provide some
    examples:

    .. code-block:: python

        # pure single excitations (equivalent options):
        uccs = UCC(excitations='s')
        uccs = UCC(excitations=1)
        uccs = UCC(excitations=[1])

        # pure double excitations (equivalent options):
        uccd = UCC(excitations='d')
        uccd = UCC(excitations=2)
        uccd = UCC(excitations=[2])

        # combinations of excitations:
        custom_ucc_sd = UCC(excitations='sd')  # see also the convenience sub-class UCCSD
        custom_ucc_sd = UCC(excitations=[1, 2])  # see also the convenience sub-class UCCSD
        custom_ucc_sdt = UCC(excitations='sdt')
        custom_ucc_sdt = UCC(excitations=[1, 2, 3])
        custom_ucc_st = UCC(excitations='st')
        custom_ucc_st = UCC(excitations=[1, 3])

        # you can even define a fully custom list of excitations:

        def custom_excitation_list(num_spin_orbitals: int,
                                   num_particles: Tuple[int, int]
                                   ) -> List[Tuple[Tuple[Any, ...], ...]]:
            # generate your list of excitations...
            my_excitation_list = [...]
            # For more information about the required format of the return statement, please take a
            # look at the documentation of
            # `qiskit_nature.circuit.library.ansatzes.utils.fermionic_excitation_generator`
            return my_excitation_list

        my_custom_ucc = UCC(excitations=custom_excitation_list)

    Keep in mind, that in all of the examples above we have not set any of the following keyword
    arguments, which must be specified before the Ansatz becomes usable:

    - `qubit_converter`
    - `num_particles`
    - `num_spin_orbitals`

    If you are using this Ansatz with a Qiskit Nature algorithm, these arguments will be set for
    you, depending on the rest of the stack.


    References:

        [1] https://arxiv.org/abs/1805.04340
    """

    EXCITATION_TYPE = {
        "s": 1,
        "d": 2,
        "t": 3,
        "q": 4,
    }

    def __init__(
        self,
        qubit_converter: Optional[QubitConverter] = None,
        num_particles: Optional[Tuple[int, int]] = None,
        num_spin_orbitals: Optional[int] = None,
        excitations: Optional[
            Union[
                str,
                int,
                List[int],
                Callable[
                    [int, Tuple[int, int]],
                    List[Tuple[Tuple[int, ...], Tuple[int, ...]]],
                ],
            ]
        ] = None,
        alpha_spin: bool = True,
        beta_spin: bool = True,
        max_spin_excitation: Optional[int] = None,
        generalized: bool = False,
        preserve_spin: bool = True,
        reps: int = 1,
        initial_state: Optional[QuantumCircuit] = None,
        evolution: Optional[EvolutionBase] = None,
    ):
        """

        Args:
            qubit_converter: the QubitConverter instance which takes care of mapping a
                :class:`~.SecondQuantizedOp` to a :class:`PauliSumOp` as well as performing all
                configured symmetry reductions on it.
            num_particles: the tuple of the number of alpha- and beta-spin particles.
            num_spin_orbitals: the number of spin orbitals.
            excitations: this can be any of the following types:

                :`str`: which contains the types of excitations. Allowed characters are
                    + `s` for singles
                    + `d` for doubles
                    + `t` for triples
                    + `q` for quadruples
                :`int`: a single, positive integer which denotes the number of excitations
                    (1 == `s`, etc.)
                :`List[int]`: a list of positive integers generalizing the above
                :`Callable`: a function which is used to generate the excitations.
                    The callable must take the __keyword__ arguments `num_spin_orbitals` and
                    `num_particles` (with identical types to those explained above) and must return
                    a `List[Tuple[Tuple[int, ...], Tuple[int, ...]]]`. For more information on how
                    to write such a callable refer to the default method
                    :meth:`~qiskit_nature.circuit.library.ansatzes.utils.generate_fermionic_excitations`.
            alpha_spin: boolean flag whether to include alpha-spin excitations.
            beta_spin: boolean flag whether to include beta-spin excitations.
            max_spin_excitation: the largest number of excitations within a spin. E.g. you can set
                this to 1 and `num_excitations` to 2 in order to obtain only mixed-spin double
                excitations (alpha,beta) but no pure-spin double excitations (alpha,alpha or
                beta,beta).
            generalized: boolean flag whether or not to use generalized excitations, which ignore
                the occupation of the spin orbitals. As such, the set of generalized excitations is
                only determined from the number of spin orbitals and independent from the number of
                particles.
            preserve_spin: boolean flag whether or not to preserve the particle spins.
            reps: The number of times to repeat the evolved operators.
            initial_state: A `QuantumCircuit` object to prepend to the circuit. Note that this
                setting does _not_ influence the `excitations`. When relying on the default
                generation method (i.e. not providing a `Callable` to `excitations`), these will
                always be constructed with respect to a `HartreeFock` reference state.
            evolution: An EvolutionBase to indicate what type of conversion to perform when 
                translating operators to gates
        """
        self._qubit_converter = qubit_converter
        self._num_particles = num_particles
        self._num_spin_orbitals = num_spin_orbitals
        self._excitations = excitations
        self._alpha_spin = alpha_spin
        self._beta_spin = beta_spin
        self._max_spin_excitation = max_spin_excitation
        self._generalized = generalized
        self._preserve_spin = preserve_spin

        if evolution is None:
            evolution = PauliTrotterEvolution()

        super().__init__(reps=reps, evolution=evolution, initial_state=initial_state)

        # To give read access to the excitation list that UCC is using.
        self._excitation_list: List[Tuple[Tuple[int, ...], Tuple[int, ...]]] = None

        # We cache these, because the generation may be quite expensive (depending on the generator)
        # and the user may want quick access to inspect these. Also, it speeds up testing for the
        # same reason!
        self._excitation_ops: List[SecondQuantizedOp] = None

    @property
    def qubit_converter(self) -> QubitConverter:
        """The qubit operator converter."""
        return self._qubit_converter

    @qubit_converter.setter
    def qubit_converter(self, conv: QubitConverter) -> None:
        """Sets the qubit operator converter."""
        self._operators = None
        self._invalidate()
        self._qubit_converter = conv

    @property
    def num_spin_orbitals(self) -> int:
        """The number of spin orbitals."""
        return self._num_spin_orbitals

    @num_spin_orbitals.setter
    def num_spin_orbitals(self, n: int) -> None:
        """Sets the number of spin orbitals."""
        self._operators = None
        self._invalidate()
        self._num_spin_orbitals = n

    @property
    def num_particles(self) -> Tuple[int, int]:
        """The number of particles."""
        return self._num_particles

    @num_particles.setter
    def num_particles(self, n: Tuple[int, int]) -> None:
        """Sets the number of particles."""
        self._operators = None
        self._invalidate()
        self._num_particles = n

    @property
    def excitations(self) -> Union[str, int, List[int], Callable]:
        """The excitations."""
        return self._excitations

    @excitations.setter
    def excitations(self, exc: Union[str, int, List[int], Callable]) -> None:
        """Sets the excitations."""
        self._operators = None
        self._invalidate()
        self._excitations = exc

    @property
    def excitation_list(self) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
        """The excitation list that UCC is using, in the required format.

        Raises:
            QiskitNatureError: If private excitation list is None.

        """
        if self._excitation_list is None:
            raise QiskitNatureError(
                "The excitation list is None. Build the operators to construct it."
            )
        return self._excitation_list

    @EvolvedOperatorAnsatz.operators.getter
    def operators(self):  # pylint: disable=invalid-overridden-method
        """The operators that are evolved in this circuit.

        Returns:
            list: The operators to be evolved contained in this ansatz or
                  None if the configuration is not complete
        """
        # Overriding the getter to build the operators on demand when they are
        # requested, if they are still set to None.
        operators = super(custom_UCC, self.__class__).operators.__get__(self)

        if operators is None or operators == [None]:
            # If the operators are None build them out if the ucc config checks out ok, otherwise
            # they will be left as None to be built at some later time.
            if self._check_ucc_configuration(raise_on_failure=False):
                # The qubit operators are cached by `EvolvedOperatorAnsatz` class. We only generate
                # them from the `SecondQuantizedOp`s produced by the generators, if they're not already
                # present. This behavior also enables the adaptive usage of the `UCC` class by
                # algorithms such as `AdaptVQE`.
                excitation_ops = self.excitation_ops()

                logger.debug("Converting SecondQuantizedOps into PauliSumOps...")
                # Convert operators according to saved state in converter from the conversion of the
                # main operator since these need to be compatible. If Z2 Symmetry tapering was done
                # it may be that one or more excitation operators do not commute with the
                # symmetry. Normally the converted operators are maintained at the same index by
                # the converter inserting None as the result if an operator did not commute. Here
                # we are not interested in that just getting the valid set of operators so that
                # behavior is suppressed.
                self.operators = self.qubit_converter.convert_match(
                    excitation_ops, suppress_none=True
                )

        return super(custom_UCC, self.__class__).operators.__get__(self)

    def _invalidate(self):
        self._excitation_ops = None
        super()._invalidate()

    def _check_configuration(self, raise_on_failure: bool = True) -> bool:
        # Check our local config is valid first. The super class will check the
        # operators by getting them, and if we detect they are still None they
        # will be built so that its valid check will end up passing in that regard.
        if not self._check_ucc_configuration(raise_on_failure):
            return False

        return super()._check_configuration(raise_on_failure)

    # pylint: disable=too-many-return-statements
    def _check_ucc_configuration(self, raise_on_failure: bool = True) -> bool:
        # Check the local config, separated out that it can be checked via build
        # or ahead of building operators to make sure everything needed is present.
        if self.num_spin_orbitals is None:
            if raise_on_failure:
                raise ValueError("The number of spin orbitals cannot be 'None'.")
            return False

        if self.num_spin_orbitals <= 0:
            if raise_on_failure:
                raise ValueError(
                    f"The number of spin orbitals must be > 0 was {self.num_spin_orbitals}."
                )
            return False

        if self.num_particles is None:
            if raise_on_failure:
                raise ValueError("The number of particles cannot be 'None'.")
            return False

        if any(n < 0 for n in self.num_particles):
            if raise_on_failure:
                raise ValueError(
                    f"The number of particles cannot be smaller than 0 was {self.num_particles}."
                )
            return False

        if sum(self.num_particles) >= self.num_spin_orbitals:
            if raise_on_failure:
                raise ValueError(
                    f"The number of spin orbitals {self.num_spin_orbitals}"
                    f"must be greater than total number of particles "
                    f"{sum(self.num_particles)}."
                )
            return False

        if self.excitations is None:
            if raise_on_failure:
                raise ValueError("The excitations cannot be `None`.")
            return False

        if self.qubit_converter is None:
            if raise_on_failure:
                raise ValueError("The qubit_converter cannot be `None`.")
            return False

        return True

    def excitation_ops(self) -> List[SecondQuantizedOp]:
        """Parses the excitations and generates the list of operators.

        Raises:
            QiskitNatureError: if invalid excitations are specified.

        Returns:
            The list of generated excitation operators.
        """
        if self._excitation_ops is not None:
            return self._excitation_ops

        excitation_list = self._get_excitation_list()

        self._check_excitation_list(excitation_list)

        logger.debug("Converting excitations into SecondQuantizedOps...")
        excitation_ops = self._build_fermionic_excitation_ops(excitation_list)

        self._excitation_list = excitation_list
        self._excitation_ops = excitation_ops
        return excitation_ops

    def _get_excitation_list(self) -> List[Tuple[Tuple[int, ...], Tuple[int, ...]]]:
        generators = self._get_excitation_generators()

        logger.debug("Generating excitation list...")
        excitations = []
        for gen in generators:
            excitations.extend(
                gen(
                    num_spin_orbitals=self.num_spin_orbitals,
                    num_particles=self.num_particles,
                    num_sub=[self.num_spin_orbitals//4, self.num_spin_orbitals//4]  # Warning: assuming two equal fragments
                )
            )

        return excitations

    def _get_excitation_generators(self) -> List[Callable]:
        logger.debug("Gathering excitation generators...")
        generators: List[Callable] = []

        extra_kwargs = {
            "alpha_spin": self._alpha_spin,
            "beta_spin": self._beta_spin,
            "max_spin_excitation": self._max_spin_excitation,
            "generalized": self._generalized,
            "preserve_spin": self._preserve_spin,
        }

        if isinstance(self.excitations, str):
            for exc in self.excitations:
                generators.append(
                    partial(
                        generate_fermionic_excitations,
                        num_excitations=self.EXCITATION_TYPE[exc],
                        **extra_kwargs,
                    )
                )
        elif isinstance(self.excitations, int):
            generators.append(
                partial(
                    generate_fermionic_excitations, num_excitations=self.excitations, **extra_kwargs
                )
            )
        elif isinstance(self.excitations, list):
            for exc in self.excitations:  # type: ignore
                generators.append(
                    partial(generate_fermionic_excitations, num_excitations=exc, **extra_kwargs)
                )
        elif callable(self.excitations):
            generators = [self.excitations]
        else:
            raise QiskitNatureError(f"Invalid excitation configuration: {self.excitations}")

        return generators

    def _check_excitation_list(self, excitations: Sequence) -> None:
        """Checks the format of the given excitation operators.

        The following conditions are checked:
        - the list of excitations consists of pairs of tuples
        - each pair of excitation indices has the same length
        - the indices within each excitation pair are unique

        Args:
            excitations: the list of excitations

        Raises:
            QiskitNatureError: if format of excitations is invalid
        """
        logger.debug("Checking excitation list...")

        error_message = "{error} in the following UCC excitation: {excitation}"

        for excitation in excitations:
            if len(excitation) != 2:
                raise QiskitNatureError(
                    error_message.format(error="Invalid number of tuples", excitation=excitation)
                    + "; Two tuples are expected, e.g. ((0, 1, 4), (2, 3, 6))"
                )

            if len(excitation[0]) != len(excitation[1]):
                raise QiskitNatureError(
                    error_message.format(
                        error="Different number of occupied and virtual indices",
                        excitation=excitation,
                    )
                )

            # Commenting this out to be able to match Matt's amplitudes
            #if any(i in excitation[0] for i in excitation[1]) or any(
            #    len(set(indices)) != len(indices) for indices in excitation
            #):
            #    raise QiskitNatureError(
            #        error_message.format(error="Duplicated indices", excitation=excitation)
            #    )

    def _build_fermionic_excitation_ops(self, excitations: Sequence) -> List[FermionicOp]:
        """Builds all possible excitation operators with the given number of excitations for the
        specified number of particles distributed in the number of orbitals.

        Args:
            excitations: the list of excitations.

        Returns:
            The list of excitation operators in the second quantized formalism.
        """
        operators = []

        for exc in excitations:
            label = ["I"] * self.num_spin_orbitals
            for occ in exc[0]:
                label[occ] = "+"
            for unocc in exc[1]:
                label[unocc] = "-"
            op = FermionicOp("".join(label), display_format="dense")
            op -= op.adjoint()
            # we need to account for an additional imaginary phase in the exponent (see also
            # `PauliTrotterEvolution.convert`)
            op *= 1j  # type: ignore
            operators.append(op)

        return operators
