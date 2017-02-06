Program shuffle

! Program to deal 4 hands as if for a game of bridge
! also sort the hands
     
  IMPLICIT NONE

  Integer, Dimension( 1:13, 1:4 ) :: hands
  Integer :: i
  TYPE card
     CHARACTER(len=1) :: suit, value
  END TYPE card
  TYPE bridge_hand
     TYPE(card), dimension(13) :: north, south, east, west
  END TYPE bridge_hand


! Generate 4 unsorted hands
  Call shuffle_cards( hands )

! Sort and write out each hand in turn
  Do i = 1, 4
     Call sort( hands( :, i ) )
     Call write_hand( i, hands( :, i ) )
  End Do
  
Contains

 
  Subroutine shuffle_cards( hands )

! Shuffle a deck of 52 cards and generate 4 hands.
! The method is to pick a card at random and store it.
! Then that card is moved to the high end of the list
! So when 2 cards have been used cards 1-50 will be the
! ones we can choose from, cards 51, 52 will be the
! ones we have already used.

    Integer, Dimension( :, : ), Intent( Out ) :: hands

    Integer :: i, card_selected
    Integer :: this_hand, card_in_hand
    Integer, Dimension(52) :: cards
    Integer, Dimension(52) :: tmp
    Real :: a_number

! Give each card in the deck a unique label
    cards = (/ (i,i=1,52) /)

! Which hand we are currently dealing to
    this_hand = 0
! Which card in the current hand we are dealing
    card_in_hand = 0
    Do i = 52,1,-1

! Every 14th card move onto a new hand
       If(Mod(i,13).Eq.0)Then
          this_hand = this_hand + 1
          card_in_hand = 0
       End If
       
! Pick a card at random and store it
       Call Random_number(a_number)
       card_selected = Int(a_number*Real(i))+1
       card_in_hand = card_in_hand + 1
       hands( card_in_hand, this_hand ) = cards( card_selected )

! Move the card into the list of cards already used
       If(card_selected.Ne.i)Then
          cards(card_selected) = cards(i)
       End If

       tmp( i ) = cards( card_selected ) 

    End Do

! Store the last card
    hands( card_in_hand, this_hand ) = cards( 1 )

  End Subroutine shuffle_cards

  Subroutine sort( a )

! Integer sorting routine. Method is selection sort.

    Integer, Dimension( : ), Intent( InOut ) :: a

    Integer, Dimension( 1:1 ) :: location

    Integer :: n
    Integer :: swap
    Integer :: i, j

    n = Size( a )

    Do i = 1, n - 1
!
! Note two things; firstly minloc returns an array.
! Secondly in returns the location of the minimum within the
! array segment passed to it, not within the whole array. Hence
! we need to add (i-1) to its value to take account of this offset
!
       location = Minloc( a( i:n ) )
       j = location( 1 ) + i - 1
       swap   = a( j )
       a( j ) = a( i )
       a( i ) = swap

    End Do

  End Subroutine sort

  Subroutine write_hand( hand_number, hand )

! Write out a hand of cards
    
    Integer                , Intent( In ) :: hand_number
    Integer, Dimension( : ), Intent( In ) :: hand
    
    Integer :: i, selected_value, selected_suit
    
    Character(len=1),Dimension(4)  :: suit
    Character(len=1),Dimension(13) :: value
    
    suit = (/ 'C','D','H','S' /)
    value = (/ '2','3','4','5','6','7','8','9','T','J','Q','K', 'A' /)
    
    Write(*,*)'hand', hand_number
    Write(*,*)'***************'
    
    Do i = 1, Size( hand )
       
       selected_suit = (hand(i)-1)/13 + 1
       selected_value = Mod(hand(i)-1,13) + 1
       Write(*,*)value(selected_value),suit(selected_suit)
       
    End Do
    
  End Subroutine write_hand
  
End Program shuffle
    
