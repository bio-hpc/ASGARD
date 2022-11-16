$mol = bless( {
  'atoms' => [
    bless( {
      'Z' => 6,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '-1.0135',
            '-0.2651',
            '-0.2957'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a1',
      'parent' => {},
      'symbol' => 'C'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 9,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '-1.7202',
            '-0.2677',
            '0.8572'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a2',
      'parent' => {},
      'symbol' => 'F'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 17,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '-1.4656',
            '1.1797',
            '-1.2790'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a3',
      'parent' => {},
      'symbol' => 'Cl'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 6,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '0.4955',
            '-0.1819',
            '-0.0726'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a4',
      'parent' => {},
      'symbol' => 'C'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 8,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '1.3585',
            '-0.6977',
            '-0.7542'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a5',
      'parent' => {},
      'symbol' => 'O'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 8,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '0.8778',
            '0.5627',
            '0.9824'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a6',
      'parent' => {},
      'symbol' => 'O'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 1,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '-1.2888',
            '-1.1841',
            '-0.8559'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a7',
      'parent' => {},
      'symbol' => 'H'
    }, 'Chemistry::Atom' ),
    bless( {
      'Z' => 1,
      'bonds' => [],
      'coords' => bless( [
        [
          [
            '1.8276',
            '0.5809',
            '1.0507'
          ]
        ],
        1,
        3
      ], 'Math::VectorReal' ),
      'id' => 'a8',
      'parent' => {},
      'symbol' => 'H'
    }, 'Chemistry::Atom' )
  ],
  'bonds' => [],
  'byId' => {
    'a1' => {},
    'a2' => {},
    'a3' => {},
    'a4' => {},
    'a5' => {},
    'a6' => {},
    'a7' => {},
    'a8' => {}
  },
  'id' => 'mol1',
  'name' => 'build.mol'
}, 'Chemistry::Mol' );
$mol->{'atoms'}[0]{'parent'} = $mol;
$mol->{'atoms'}[1]{'parent'} = $mol;
$mol->{'atoms'}[2]{'parent'} = $mol;
$mol->{'atoms'}[3]{'parent'} = $mol;
$mol->{'atoms'}[4]{'parent'} = $mol;
$mol->{'atoms'}[5]{'parent'} = $mol;
$mol->{'atoms'}[6]{'parent'} = $mol;
$mol->{'atoms'}[7]{'parent'} = $mol;
$mol->{'byId'}{'a1'} = $mol->{'atoms'}[0];
$mol->{'byId'}{'a2'} = $mol->{'atoms'}[1];
$mol->{'byId'}{'a3'} = $mol->{'atoms'}[2];
$mol->{'byId'}{'a4'} = $mol->{'atoms'}[3];
$mol->{'byId'}{'a5'} = $mol->{'atoms'}[4];
$mol->{'byId'}{'a6'} = $mol->{'atoms'}[5];
$mol->{'byId'}{'a7'} = $mol->{'atoms'}[6];
$mol->{'byId'}{'a8'} = $mol->{'atoms'}[7];
